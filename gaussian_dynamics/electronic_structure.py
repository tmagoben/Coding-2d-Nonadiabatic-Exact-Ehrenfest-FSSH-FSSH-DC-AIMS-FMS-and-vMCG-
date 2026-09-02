from dataclasses import dataclass, field
from typing import Protocol
import hashlib
import json
import numpy as np

from .adiabatic import adiabatic_point


@dataclass
class ElectronicStructurePoint:
    """Electronic data resolved along one generalized nuclear coordinate q."""
    q: float
    energies: np.ndarray
    gradients_q: np.ndarray
    nac_q: np.ndarray
    metadata: dict = field(default_factory=dict)

    def validate(self):
        e=np.asarray(self.energies,float)
        g=np.asarray(self.gradients_q,float)
        d=np.asarray(self.nac_q,float)

        if e.ndim != 1:
            raise ValueError("energies must be one-dimensional.")
        if g.shape != e.shape:
            raise ValueError("gradients_q must have one value per state.")
        if d.shape != (len(e),len(e)):
            raise ValueError("nac_q must have shape (nstates,nstates).")
        if not np.all(np.isfinite(e)) or not np.all(np.isfinite(g)) or not np.all(np.isfinite(d)):
            raise ValueError("Electronic structure data contain non-finite values.")
        if not np.allclose(d, -d.T, atol=1e-8):
            raise ValueError("For this real-state interface, nac_q must be antisymmetric.")
        return self


class ElectronicStructureProvider(Protocol):
    def evaluate(self, q: float) -> ElectronicStructurePoint:
        ...


class AnalyticAvoidedCrossingProvider:
    """Provider exposing the analytic v0.2 avoided-crossing model."""
    def evaluate(self, q):
        E,U,g,d=adiabatic_point(float(q))
        return ElectronicStructurePoint(
            q=float(q),
            energies=E.copy(),
            gradients_q=g.copy(),
            nac_q=d.copy(),
            metadata={"provider":"analytic_avoided_crossing"},
        ).validate()


class TabulatedElectronicStructureProvider:
    """Linear interpolation of prevalidated 1D electronic-structure data.

    Input data must already be state tracked and gauge consistent.
    This class deliberately does not silently reorder electronic states.
    """
    def __init__(self, q_grid, energies, gradients_q, nac_q):
        self.q=np.asarray(q_grid,float)
        self.E=np.asarray(energies,float)
        self.G=np.asarray(gradients_q,float)
        self.D=np.asarray(nac_q,float)

        if np.any(np.diff(self.q) <= 0):
            raise ValueError("q_grid must be strictly increasing.")
        n=len(self.q)
        if self.E.shape[0] != n or self.G.shape[0] != n or self.D.shape[0] != n:
            raise ValueError("All tables must use the same q_grid.")
        if self.E.shape != self.G.shape:
            raise ValueError("energies and gradients must have matching shapes.")
        ns=self.E.shape[1]
        if self.D.shape != (n,ns,ns):
            raise ValueError("nac_q has incompatible shape.")

        for i in range(n):
            ElectronicStructurePoint(self.q[i],self.E[i],self.G[i],self.D[i]).validate()

    def _interp(self, arr, q):
        q=float(q)
        if q < self.q[0] or q > self.q[-1]:
            raise ValueError("Requested geometry lies outside the tabulated domain.")
        if arr.ndim == 2:
            return np.array([np.interp(q,self.q,arr[:,j]) for j in range(arr.shape[1])])
        if arr.ndim == 3:
            ns=arr.shape[1]
            out=np.zeros((ns,ns))
            for i in range(ns):
                for j in range(ns):
                    out[i,j]=np.interp(q,self.q,arr[:,i,j])
            return out
        raise ValueError("Unsupported table rank.")

    def evaluate(self,q):
        point=ElectronicStructurePoint(
            q=float(q),
            energies=self._interp(self.E,q),
            gradients_q=self._interp(self.G,q),
            nac_q=self._interp(self.D,q),
            metadata={"provider":"tabulated_linear_interpolation"},
        )
        # Reimpose exact antisymmetry after interpolation.
        point.nac_q=0.5*(point.nac_q-point.nac_q.T)
        return point.validate()


class CachedProvider:
    """Deterministic in-memory cache around any provider."""
    def __init__(self, provider, digits=12):
        self.provider=provider
        self.digits=int(digits)
        self.cache={}
        self.hits=0
        self.misses=0

    def _key(self,q):
        return round(float(q),self.digits)

    def evaluate(self,q):
        key=self._key(q)
        if key in self.cache:
            self.hits += 1
            return self.cache[key]
        self.misses += 1
        point=self.provider.evaluate(key)
        self.cache[key]=point
        return point


def project_cartesian_vector_to_coordinate(vector_R, tangent_dR_dq):
    """Chain-rule projection A_q = sum_A A_R_A dot dR_A/dq."""
    v=np.asarray(vector_R,float)
    t=np.asarray(tangent_dR_dq,float)
    if v.shape != t.shape:
        raise ValueError("Cartesian vector and coordinate tangent must have equal shape.")
    return float(np.sum(v*t))


def point_fingerprint(point):
    """SHA-256 of the numerical point + metadata for provenance/regression use."""
    payload={
        "q":float(point.q),
        "energies":np.asarray(point.energies).tolist(),
        "gradients_q":np.asarray(point.gradients_q).tolist(),
        "nac_q":np.asarray(point.nac_q).tolist(),
        "metadata":point.metadata,
    }
    raw=json.dumps(payload,sort_keys=True,separators=(",",":")).encode()
    return hashlib.sha256(raw).hexdigest()
