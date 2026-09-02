from dataclasses import dataclass, field
import hashlib
import json
import time
import warnings
import numpy as np

from .molecular_backend import (
    GeneralizedElectronicStructurePoint,
    geometry_fingerprint,
)
from .molecular_snapshot_v19 import (
    MolecularElectronicSnapshotV19,
    TrackedGeneralizedSnapshotV19,
)
from .state_tracking import (
    transform_state_properties,
)
from .state_tracking_v19 import (
    scalable_maximum_overlap_assignment_v19,
)


@dataclass(frozen=True)
class BackendEvaluationPolicyV19:
    retries: int=0
    failure_policy: str="raise"  # raise | nearest_cache
    max_fallback_distance: float=0.05

    def validate(self):
        if self.retries<0:
            raise ValueError("retries cannot be negative.")
        if self.failure_policy not in {"raise","nearest_cache"}:
            raise ValueError(
                "failure_policy must be 'raise' or 'nearest_cache'."
            )
        if self.max_fallback_distance<0.0:
            raise ValueError("max_fallback_distance cannot be negative.")
        return self


@dataclass(frozen=True)
class MolecularTrackingSettingsV19:
    minimum_overlap: float=0.50
    minimum_score_margin: float=0.05
    ambiguity_policy: str="raise"  # raise | warn | accept
    real_gauge: bool=True
    cache_digits: int=12
    reference_policy: str="nearest"  # currently nearest only

    def validate(self):
        if self.ambiguity_policy not in {"raise","warn","accept"}:
            raise ValueError("invalid ambiguity_policy.")
        if self.reference_policy!="nearest":
            raise ValueError("v0.19 currently supports reference_policy='nearest'.")
        if self.cache_digits<0:
            raise ValueError("cache_digits cannot be negative.")
        return self


@dataclass
class MolecularProviderDiagnosticsV19:
    evaluate_calls: int=0
    cache_hits: int=0
    cache_misses: int=0
    backend_attempts: int=0
    backend_failures: int=0
    fallback_uses: int=0
    tracking_ambiguities: int=0
    total_backend_seconds: float=0.0
    maximum_backend_seconds: float=0.0
    history: list=field(default_factory=list)

    def as_dict(self):
        return {
            "evaluate_calls":int(self.evaluate_calls),
            "cache_hits":int(self.cache_hits),
            "cache_misses":int(self.cache_misses),
            "backend_attempts":int(self.backend_attempts),
            "backend_failures":int(self.backend_failures),
            "fallback_uses":int(self.fallback_uses),
            "tracking_ambiguities":int(self.tracking_ambiguities),
            "total_backend_seconds":float(self.total_backend_seconds),
            "maximum_backend_seconds":float(self.maximum_backend_seconds),
            "history":list(self.history),
        }


def _tracking_overlap(previous, current, overlap_engine=None):
    if overlap_engine is not None:
        return np.asarray(
            overlap_engine(previous,current),
            dtype=complex,
        )

    if (
        previous.state_vectors is not None
        and current.state_vectors is not None
    ):
        return (
            previous.state_vectors.conj().T
            @current.state_vectors
        )

    raise ValueError(
        "No cross-geometry overlap representation is available. "
        "Provide state_vectors or an overlap_engine."
    )


def _transform_snapshot(raw, result):
    point=raw.point
    E,G,D=transform_state_properties(
        point.energies,
        point.gradients_cart,
        point.nac_cart,
        result,
    )

    if raw.state_vectors is not None:
        perm=np.asarray(result.permutation,dtype=int)
        phase=np.asarray(result.phase_factors,dtype=complex)
        vectors=raw.state_vectors[:,perm]*phase[None,:]
    else:
        vectors=None

    wave=raw.wavefunction_snapshot
    if wave is not None and hasattr(wave,"with_transformed_roots"):
        wave=wave.with_transformed_roots(
            result.permutation,
            result.phase_factors,
        )

    from .molecular_backend import CartesianElectronicStructurePoint
    tracked_point=CartesianElectronicStructurePoint(
        geometry=point.geometry,
        energies=np.asarray(E,float),
        gradients_cart=np.asarray(G,float),
        nac_cart=np.asarray(D,float),
        masses_amu=np.asarray(point.masses_amu,float),
        scaled_nac_cart=(
            None
            if point.scaled_nac_cart is None
            else np.asarray(point.scaled_nac_cart).copy()
        ),
        metadata={
            **dict(point.metadata),
            "v19_tracking":result.as_metadata(),
        },
    ).validate()

    return MolecularElectronicSnapshotV19(
        point=tracked_point,
        state_vectors=vectors,
        wavefunction_snapshot=wave,
        metadata={
            **dict(raw.metadata),
            "tracking":result.as_metadata(),
        },
    ).validate()


class TrackedMolecularDirectProviderV19:
    """Order-tolerant cached molecular provider using nearest-anchor overlap tracking.

    Unlike the sequential v0.6 tracked backend, each new generalized-coordinate point
    is aligned to the nearest already accepted cached node. This is better suited to
    branched TBF centers and pair centroids, where evaluation order is not a physical
    trajectory ordering.

    It remains a local discrete tracking rule, not a proof of globally path-independent
    gauge transport around nontrivial holonomy.
    """

    def __init__(
        self,
        backend,
        geometry_map,
        *,
        overlap_engine=None,
        tracking=MolecularTrackingSettingsV19(),
        failure=BackendEvaluationPolicyV19(),
    ):
        self.backend=backend
        self.geometry_map=geometry_map
        self.overlap_engine=overlap_engine
        self.tracking=tracking.validate()
        self.failure=failure.validate()
        self._cache={}
        self._q={}
        self._snapshots={}
        self.diagnostics=MolecularProviderDiagnosticsV19()

    def _key(self,q):
        q=np.asarray(q,float)
        return tuple(
            round(float(x),self.tracking.cache_digits)
            for x in q
        )

    def _node_id(self,key):
        raw=json.dumps(key,separators=(",",":")).encode()
        return "mol-"+hashlib.sha256(raw).hexdigest()[:16]

    def _after_cache_insert(self,key,q):
        """Extension hook for indexed cache implementations."""
        return None

    def _nearest_key(self,q):
        if not self._q:
            return None,None
        q=np.asarray(q,float)
        keys=list(self._q)
        distances=np.asarray([
            np.linalg.norm(q-self._q[k])
            for k in keys
        ])
        m=float(np.min(distances))
        candidates=[
            keys[i] for i,d in enumerate(distances)
            if abs(float(d)-m)<=1e-14
        ]
        key=sorted(candidates,key=repr)[0]
        return key,m

    def _backend_snapshot(self,geometry):
        last=None
        for _ in range(self.failure.retries+1):
            self.diagnostics.backend_attempts+=1
            t0=time.perf_counter()
            try:
                if hasattr(self.backend,"evaluate_snapshot"):
                    snap=self.backend.evaluate_snapshot(geometry)
                else:
                    point=self.backend.evaluate(geometry)
                    raise TypeError(
                        "Backend lacks evaluate_snapshot(); an overlap-capable "
                        "snapshot adapter is required for v0.19 tracking."
                    )
                elapsed=time.perf_counter()-t0
                self.diagnostics.total_backend_seconds+=elapsed
                self.diagnostics.maximum_backend_seconds=max(
                    self.diagnostics.maximum_backend_seconds,
                    elapsed,
                )
                return snap.validate(),elapsed
            except Exception as exc:
                elapsed=time.perf_counter()-t0
                self.diagnostics.total_backend_seconds+=elapsed
                self.diagnostics.maximum_backend_seconds=max(
                    self.diagnostics.maximum_backend_seconds,
                    elapsed,
                )
                self.diagnostics.backend_failures+=1
                last=exc
        raise last

    def _project(self,q,snapshot):
        point=snapshot.point
        J=self.geometry_map.J
        ns=len(point.energies)
        grad_flat=np.asarray(
            point.gradients_cart,float
        ).reshape(ns,-1)
        nac_flat=np.asarray(
            point.nac_cart,float
        ).reshape(ns,ns,-1)

        gradients_q=grad_flat@J
        nac_q=np.einsum(
            "ijr,ra->ija",
            nac_flat,J,
        )
        Mq=self.geometry_map.mass_matrix_q_au(
            point.masses_amu
        )
        return GeneralizedElectronicStructurePoint(
            q=np.asarray(q,float).copy(),
            energies=np.asarray(point.energies,float).copy(),
            gradients_q=gradients_q,
            nac_q=nac_q,
            mass_matrix_q_au=Mq,
            metadata=dict(point.metadata),
        ).validate()

    def _fallback(self,q,failed_exc):
        nearest,distance=self._nearest_key(q)
        if (
            nearest is None
            or distance>self.failure.max_fallback_distance
        ):
            raise failed_exc

        source=self._cache[nearest]
        p=source.point
        fallback=GeneralizedElectronicStructurePoint(
            q=np.asarray(q,float).copy(),
            energies=p.energies.copy(),
            gradients_q=p.gradients_q.copy(),
            nac_q=p.nac_q.copy(),
            mass_matrix_q_au=p.mass_matrix_q_au.copy(),
            metadata={
                **dict(p.metadata),
                "v19_failure_fallback":True,
                "fallback_source_node":
                    source.node_id,
                "fallback_distance":float(distance),
                "fallback_exception":
                    type(failed_exc).__name__,
            },
        ).validate()
        self.diagnostics.fallback_uses+=1
        self.diagnostics.history.append({
            "kind":"fallback",
            "q":np.asarray(q,float).tolist(),
            "source":source.node_id,
            "distance":float(distance),
        })
        return fallback

    def evaluate_snapshot(self,q):
        q=np.asarray(q,float)
        if q.shape!=(self.geometry_map.nq,):
            raise ValueError("q has incompatible shape.")
        self.diagnostics.evaluate_calls+=1
        key=self._key(q)

        if key in self._cache:
            self.diagnostics.cache_hits+=1
            self.diagnostics.history.append({
                "kind":"cache_hit",
                "q":q.tolist(),
                "node":self._cache[key].node_id,
            })
            return self._cache[key]

        self.diagnostics.cache_misses+=1
        geometry=self.geometry_map.geometry(q)

        try:
            raw,elapsed=self._backend_snapshot(geometry)
        except Exception as exc:
            if self.failure.failure_policy=="nearest_cache":
                point=self._fallback(q,exc)
                node=self._node_id(key)
                return TrackedGeneralizedSnapshotV19(
                    point=point,
                    state_vectors=None,
                    wavefunction_snapshot=None,
                    node_id=node,
                    tracking_metadata={
                        "fallback":True,
                    },
                    source_metadata={},
                )
            raise

        reference_key,distance=self._nearest_key(q)

        if reference_key is None:
            tracked=raw
            tracking_meta={
                "reference_node":None,
                "reference_distance":None,
                "initial_reference":True,
                "ambiguous":False,
            }
        else:
            reference=self._snapshots[reference_key]
            overlap=_tracking_overlap(
                reference,raw,
                self.overlap_engine,
            )
            result=scalable_maximum_overlap_assignment_v19(
                overlap,
                minimum_overlap=
                    self.tracking.minimum_overlap,
                minimum_score_margin=
                    self.tracking.minimum_score_margin,
                real_gauge=
                    self.tracking.real_gauge,
            )
            if result.ambiguous:
                self.diagnostics.tracking_ambiguities+=1
                msg=(
                    "Ambiguous v0.19 molecular state tracking: "
                    +"; ".join(result.reasons)
                )
                if self.tracking.ambiguity_policy=="raise":
                    raise RuntimeError(msg)
                if self.tracking.ambiguity_policy=="warn":
                    warnings.warn(msg,RuntimeWarning)

            tracked=_transform_snapshot(raw,result)
            tracking_meta={
                **result.as_metadata(),
                "reference_node":
                    self._cache[reference_key].node_id,
                "reference_distance":
                    float(distance),
                "initial_reference":False,
            }

        point=self._project(q,tracked)
        point.metadata.update({
            "v19_node_id":self._node_id(key),
            "v19_backend_seconds":float(elapsed),
            "v19_cache_hit":False,
            "v19_tracking_reference":
                tracking_meta.get("reference_node"),
        })

        accepted=TrackedGeneralizedSnapshotV19(
            point=point,
            state_vectors=(
                None
                if tracked.state_vectors is None
                else tracked.state_vectors.copy()
            ),
            wavefunction_snapshot=
                tracked.wavefunction_snapshot,
            node_id=self._node_id(key),
            tracking_metadata=tracking_meta,
            source_metadata=dict(tracked.metadata),
        ).validate()

        self._cache[key]=accepted
        self._q[key]=q.copy()
        self._snapshots[key]=tracked
        self._after_cache_insert(key,q)

        self.diagnostics.history.append({
            "kind":"backend_miss",
            "q":q.tolist(),
            "node":accepted.node_id,
            "reference":
                tracking_meta.get("reference_node"),
            "reference_distance":
                tracking_meta.get("reference_distance"),
            "backend_seconds":float(elapsed),
        })
        return accepted

    def evaluate(self,q):
        return self.evaluate_snapshot(q).point

    def snapshot_overlap(self,left,right):
        """Return <Phi(left)|Phi(right)> in the accepted tracked frames."""
        if (
            left.state_vectors is not None
            and right.state_vectors is not None
        ):
            return (
                left.state_vectors.conj().T
                @right.state_vectors
            )
        if (
            self.overlap_engine is not None
            and left.wavefunction_snapshot is not None
            and right.wavefunction_snapshot is not None
        ):
            return np.asarray(
                self.overlap_engine(left,right),
                dtype=complex,
            )
        raise ValueError(
            "No accepted cross-geometry overlap representation is available."
        )

    def cost_estimate(
        self,
        q,
        *,
        cached_cost=0.1,
        nearby_cost=0.5,
        new_cost=5.0,
        nearby_radius=0.05,
    ):
        q=np.asarray(q,float)
        key=self._key(q)
        if key in self._cache:
            return {
                "normalized_cost":float(cached_cost),
                "cache_hit":True,
                "nearby_cache":True,
                "distance":0.0,
            }

        _,distance=self._nearest_key(q)
        if distance is not None and distance<=nearby_radius:
            return {
                "normalized_cost":float(nearby_cost),
                "cache_hit":False,
                "nearby_cache":True,
                "distance":float(distance),
            }
        return {
            "normalized_cost":float(new_cost),
            "cache_hit":False,
            "nearby_cache":False,
            "distance":None if distance is None else float(distance),
        }

    def cache_size(self):
        return len(self._cache)

    def cached_nodes(self):
        return tuple(
            self._cache[k] for k in sorted(self._cache,key=repr)
        )

    def diagnostics_dict(self):
        out=self.diagnostics.as_dict()
        out["cache_size"]=self.cache_size()
        return out
