from dataclasses import dataclass, field
from typing import Protocol
import hashlib
import json
import numpy as np


# 2022 CODATA-scale conversion used by the dynamics layer.
AMU_TO_ELECTRON_MASS = 1822.888486209


@dataclass(frozen=True)
class MolecularGeometry:
    symbols: tuple
    coords_bohr: np.ndarray

    def __post_init__(self):
        coords = np.asarray(self.coords_bohr, dtype=float)
        if coords.ndim != 2 or coords.shape[1] != 3:
            raise ValueError("coords_bohr must have shape (natom,3).")
        if len(self.symbols) != coords.shape[0]:
            raise ValueError("symbols and coords_bohr have inconsistent atom counts.")
        if not np.all(np.isfinite(coords)):
            raise ValueError("Geometry contains non-finite coordinates.")
        object.__setattr__(self, "symbols", tuple(str(s) for s in self.symbols))
        object.__setattr__(self, "coords_bohr", coords.copy())

    @property
    def natom(self):
        return len(self.symbols)


@dataclass
class CartesianElectronicStructurePoint:
    geometry: MolecularGeometry
    energies: np.ndarray
    gradients_cart: np.ndarray
    nac_cart: np.ndarray
    masses_amu: np.ndarray
    scaled_nac_cart: np.ndarray | None = None
    metadata: dict = field(default_factory=dict)

    def validate(self, atol=1e-8):
        e = np.asarray(self.energies, dtype=float)
        g = np.asarray(self.gradients_cart, dtype=float)
        d = np.asarray(self.nac_cart, dtype=float)
        m = np.asarray(self.masses_amu, dtype=float)

        ns = len(e)
        na = self.geometry.natom

        if e.ndim != 1:
            raise ValueError("energies must be one-dimensional.")
        if g.shape != (ns, na, 3):
            raise ValueError("gradients_cart must have shape (nstate,natom,3).")
        if d.shape != (ns, ns, na, 3):
            raise ValueError("nac_cart must have shape (nstate,nstate,natom,3).")
        if m.shape != (na,):
            raise ValueError("masses_amu must have shape (natom,).")
        if np.any(m <= 0.0):
            raise ValueError("Atomic masses must be positive.")

        for arr, name in [(e, "energies"), (g, "gradients"), (d, "NACs"), (m, "masses")]:
            if not np.all(np.isfinite(arr)):
                raise ValueError(f"{name} contain non-finite values.")

        if not np.allclose(d, -np.swapaxes(d, 0, 1), atol=atol):
            raise ValueError("Real nac_cart must be antisymmetric in state indices.")

        if not np.allclose(np.diagonal(d, axis1=0, axis2=1), 0.0, atol=atol):
            raise ValueError("Diagonal NACs must be zero in this real-state contract.")

        if self.scaled_nac_cart is not None:
            scaled = np.asarray(self.scaled_nac_cart, dtype=float)
            if scaled.shape != d.shape:
                raise ValueError("scaled_nac_cart has incompatible shape.")
            if not np.all(np.isfinite(scaled)):
                raise ValueError("scaled_nac_cart contains non-finite values.")
            self.scaled_nac_cart = scaled

        self.energies = e
        self.gradients_cart = g
        self.nac_cart = d
        self.masses_amu = m
        return self


class MolecularElectronicStructureBackend(Protocol):
    def evaluate(self, geometry: MolecularGeometry) -> CartesianElectronicStructurePoint:
        ...


@dataclass
class GeneralizedElectronicStructurePoint:
    q: np.ndarray
    energies: np.ndarray
    gradients_q: np.ndarray
    nac_q: np.ndarray
    mass_matrix_q_au: np.ndarray
    metadata: dict = field(default_factory=dict)

    def validate(self, atol=1e-8):
        q = np.asarray(self.q, dtype=float)
        e = np.asarray(self.energies, dtype=float)
        g = np.asarray(self.gradients_q, dtype=float)
        d = np.asarray(self.nac_q, dtype=float)
        M = np.asarray(self.mass_matrix_q_au, dtype=float)

        ns = len(e)
        nq = len(q)

        if g.shape != (ns, nq):
            raise ValueError("gradients_q must have shape (nstate,nq).")
        if d.shape != (ns, ns, nq):
            raise ValueError("nac_q must have shape (nstate,nstate,nq).")
        if M.shape != (nq, nq):
            raise ValueError("mass_matrix_q_au must have shape (nq,nq).")
        if not np.allclose(M, M.T, atol=1e-12):
            raise ValueError("mass_matrix_q_au must be symmetric.")
        if np.min(np.linalg.eigvalsh(M)) <= 0.0:
            raise ValueError("mass_matrix_q_au must be positive definite.")
        if not np.allclose(d, -np.swapaxes(d, 0, 1), atol=atol):
            raise ValueError("nac_q must be antisymmetric.")

        for arr in (q, e, g, d, M):
            if not np.all(np.isfinite(arr)):
                raise ValueError("Generalized electronic point contains non-finite data.")

        self.q = q
        self.energies = e
        self.gradients_q = g
        self.nac_q = d
        self.mass_matrix_q_au = M
        return self


class LinearGeometryMap:
    """R(q) = R0 + sum_alpha q_alpha mode_alpha."""

    def __init__(self, symbols, reference_bohr, modes):
        self.symbols = tuple(symbols)
        self.reference_bohr = np.asarray(reference_bohr, dtype=float)
        self.modes = np.asarray(modes, dtype=float)

        if self.reference_bohr.shape != (len(self.symbols), 3):
            raise ValueError("reference_bohr must have shape (natom,3).")
        if self.modes.ndim != 3 or self.modes.shape[1:] != self.reference_bohr.shape:
            raise ValueError("modes must have shape (nq,natom,3).")

        # J maps generalized coordinate displacement to flattened Cartesian displacement.
        self.J = self.modes.reshape(self.modes.shape[0], -1).T

        if np.linalg.matrix_rank(self.J) < self.J.shape[1]:
            raise ValueError("Generalized-coordinate modes are linearly dependent.")

    @property
    def nq(self):
        return self.modes.shape[0]

    def geometry(self, q):
        q = np.asarray(q, dtype=float)
        if q.shape != (self.nq,):
            raise ValueError("q has incompatible shape.")
        coords = self.reference_bohr + np.tensordot(q, self.modes, axes=(0, 0))
        return MolecularGeometry(self.symbols, coords)

    def mass_matrix_q_au(self, masses_amu):
        masses_amu = np.asarray(masses_amu, dtype=float)
        if masses_amu.shape != (len(self.symbols),):
            raise ValueError("masses_amu has incompatible shape.")
        masses_cart = np.repeat(masses_amu * AMU_TO_ELECTRON_MASS, 3)
        M_cart = np.diag(masses_cart)
        return self.J.T @ M_cart @ self.J

    @classmethod
    def cartesian(cls, geometry):
        geometry = MolecularGeometry(geometry.symbols, geometry.coords_bohr)
        ncart = 3 * geometry.natom
        modes = np.eye(ncart).reshape(ncart, geometry.natom, 3)
        return cls(geometry.symbols, geometry.coords_bohr, modes)


class GeneralizedCoordinateProvider:
    """Project a Cartesian backend onto a constant linear generalized-coordinate map."""

    def __init__(self, backend, geometry_map):
        self.backend = backend
        self.geometry_map = geometry_map

    def evaluate(self, q):
        q = np.asarray(q, dtype=float)
        geometry = self.geometry_map.geometry(q)
        point = self.backend.evaluate(geometry).validate()

        ns = len(point.energies)
        nq = self.geometry_map.nq
        J = self.geometry_map.J

        grad_flat = point.gradients_cart.reshape(ns, -1)
        nac_flat = point.nac_cart.reshape(ns, ns, -1)

        gradients_q = grad_flat @ J
        nac_q = np.einsum("ijr,ra->ija", nac_flat, J)

        Mq = self.geometry_map.mass_matrix_q_au(point.masses_amu)

        metadata = dict(point.metadata)
        metadata.update({
            "coordinate_map": "linear",
            "n_generalized_coordinates": nq,
        })

        return GeneralizedElectronicStructurePoint(
            q=q.copy(),
            energies=point.energies.copy(),
            gradients_q=gradients_q,
            nac_q=nac_q,
            mass_matrix_q_au=Mq,
            metadata=metadata,
        ).validate()


def geometry_fingerprint(geometry):
    geometry = MolecularGeometry(geometry.symbols, geometry.coords_bohr)
    payload = {
        "symbols": list(geometry.symbols),
        "coords_bohr": np.asarray(geometry.coords_bohr).round(14).tolist(),
    }
    raw = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(raw).hexdigest()


def point_fingerprint(point):
    payload = {
        "geometry": geometry_fingerprint(point.geometry),
        "energies": np.asarray(point.energies).tolist(),
        "gradients_cart": np.asarray(point.gradients_cart).tolist(),
        "nac_cart": np.asarray(point.nac_cart).tolist(),
        "masses_amu": np.asarray(point.masses_amu).tolist(),
        "metadata": point.metadata,
    }
    raw = json.dumps(payload, sort_keys=True, separators=(",", ":"), default=str).encode()
    return hashlib.sha256(raw).hexdigest()
