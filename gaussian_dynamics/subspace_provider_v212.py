from dataclasses import dataclass, field
import numpy as np

from .indexed_molecular_provider_v20 import BufferedKDTreeIndexV20
from .subspace_tracking_v21 import procrustes_subspace_alignment_v21


@dataclass(frozen=True)
class SubspaceTrackingSettingsV212:
    minimum_singular_value: float=0.85
    ambiguity_policy: str="raise"  # raise | warn | accept
    cache_digits: int=12
    rebuild_batch: int=32

    def validate(self):
        if not (0.0<=self.minimum_singular_value<=1.0):
            raise ValueError("minimum_singular_value must lie in [0,1].")
        if self.ambiguity_policy not in {"raise","warn","accept"}:
            raise ValueError("invalid ambiguity_policy.")
        return self


@dataclass
class SubspaceProviderDiagnosticsV212:
    evaluate_calls: int=0
    cache_hits: int=0
    new_points: int=0
    subspace_checks: int=0
    subspace_ambiguities: int=0
    minimum_seen_singular_value: float=1.0
    history: list=field(default_factory=list)

    def as_dict(self):
        return {
            "evaluate_calls":int(self.evaluate_calls),
            "cache_hits":int(self.cache_hits),
            "new_points":int(self.new_points),
            "subspace_checks":int(self.subspace_checks),
            "subspace_ambiguities":int(self.subspace_ambiguities),
            "minimum_seen_singular_value":float(self.minimum_seen_singular_value),
            "history":list(self.history),
        }


class SubspaceAwareOperatorProviderV212:
    """Nearest-anchor full-subspace diagnostics without forcing a gauge rotation.

    The block framework is already covariant under arbitrary local electronic frames,
    so it does not need root-by-root phase smoothing.  This provider instead checks the
    *entire represented electronic subspace* using singular values of the cross-geometry
    overlap.  It records the Procrustes transform as a diagnostic but returns the raw
    local frame unchanged.  This avoids inventing an uncomputed derivative of a
    geometry-dependent gauge transformation and preserves possible nontrivial holonomy.
    """

    def __init__(self,base_provider,dimension,settings=SubspaceTrackingSettingsV212()):
        self.base_provider=base_provider
        self.settings=settings.validate()
        self._cache={}
        self._q={}
        self.index=BufferedKDTreeIndexV20(int(dimension),self.settings.rebuild_batch)
        self.diagnostics=SubspaceProviderDiagnosticsV212()

    def _key(self,q):
        return tuple(round(float(x),self.settings.cache_digits) for x in np.asarray(q,float))

    def evaluate_snapshot(self,q):
        import warnings
        q=np.asarray(q,float)
        self.diagnostics.evaluate_calls+=1
        key=self._key(q)
        if key in self._cache:
            self.diagnostics.cache_hits+=1
            return self._cache[key]

        snap=self.base_provider.evaluate_snapshot(q)
        nearest,distance=self.index.nearest(q)
        record={"q":q.tolist(),"reference":None,"distance":None,"initial":nearest is None}
        if nearest is not None:
            ref=self._cache[nearest]
            overlap=np.asarray(self.base_provider.snapshot_overlap(ref,snap),dtype=complex)
            result=procrustes_subspace_alignment_v21(overlap)
            self.diagnostics.subspace_checks+=1
            self.diagnostics.minimum_seen_singular_value=min(
                self.diagnostics.minimum_seen_singular_value,result.minimum_singular_value
            )
            ambiguous=result.minimum_singular_value<self.settings.minimum_singular_value
            if ambiguous:
                self.diagnostics.subspace_ambiguities+=1
                msg=(
                    "Electronic subspace overlap lost continuity: minimum singular value "
                    f"{result.minimum_singular_value:.6g} < {self.settings.minimum_singular_value:.6g}."
                )
                if self.settings.ambiguity_policy=="raise":
                    raise RuntimeError(msg)
                if self.settings.ambiguity_policy=="warn":
                    warnings.warn(msg,RuntimeWarning)
            record.update({
                "reference":repr(nearest),
                "distance":float(distance),
                "minimum_singular_value":float(result.minimum_singular_value),
                "principal_angle_max":float(result.principal_angle_max),
                "ambiguous":bool(ambiguous),
                "procrustes_transform":result.transform.tolist(),
            })

        self._cache[key]=snap
        self._q[key]=q.copy()
        self.index.insert(key,q)
        self.diagnostics.new_points+=1
        self.diagnostics.history.append(record)
        return snap

    def evaluate(self,q):
        return self.evaluate_snapshot(q).point

    def snapshot_overlap(self,left,right):
        return self.base_provider.snapshot_overlap(left,right)

    def diagnostics_dict(self):
        return {
            "subspace":self.diagnostics.as_dict(),
            "spatial_index":self.index.diagnostics.as_dict(),
            "base":self.base_provider.diagnostics_dict() if hasattr(self.base_provider,"diagnostics_dict") else {},
        }
