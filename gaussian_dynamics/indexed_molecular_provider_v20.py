from dataclasses import dataclass
import numpy as np
from scipy.spatial import cKDTree

from .molecular_direct_provider_v19 import (
    TrackedMolecularDirectProviderV19,
)


@dataclass
class IndexedCacheDiagnosticsV20:
    nearest_queries: int=0
    kd_queries: int=0
    buffer_distance_checks: int=0
    rebuilds: int=0
    indexed_points: int=0
    buffered_points: int=0

    def as_dict(self):
        return {
            "nearest_queries":int(self.nearest_queries),
            "kd_queries":int(self.kd_queries),
            "buffer_distance_checks":
                int(self.buffer_distance_checks),
            "rebuilds":int(self.rebuilds),
            "indexed_points":int(self.indexed_points),
            "buffered_points":int(self.buffered_points),
        }


class BufferedKDTreeIndexV20:
    """Exact nearest-neighbor index with a bounded recent-insertion buffer.

    The immutable cKDTree is rebuilt after `rebuild_batch` trusted insertions.  Points
    added since the most recent rebuild are searched exactly by a small brute-force
    buffer.  Therefore every nearest query is exact over all trusted points while the
    normal query work is approximately

        O(log N_indexed + B*nq)

    with B <= rebuild_batch.

    Rebuild work is separately counted and should remain negligible relative to an
    ab-initio electronic-structure call for practical molecular workloads.
    """

    def __init__(self,dimension,rebuild_batch=32):
        self.dimension=int(dimension)
        self.rebuild_batch=max(int(rebuild_batch),1)
        self._all_keys=[]
        self._all_points=[]
        self._indexed_count=0
        self._tree=None
        self._indexed_keys=[]
        self.diagnostics=IndexedCacheDiagnosticsV20()

    def insert(self,key,q):
        q=np.asarray(q,dtype=float)
        if q.shape!=(self.dimension,):
            raise ValueError("indexed coordinate has incompatible shape.")
        self._all_keys.append(key)
        self._all_points.append(q.copy())
        if (
            self._tree is None
            or len(self._all_keys)-self._indexed_count
            >=self.rebuild_batch
        ):
            self.rebuild()
        else:
            self._refresh_counts()

    def rebuild(self):
        if not self._all_points:
            self._tree=None
            self._indexed_keys=[]
            self._indexed_count=0
            self._refresh_counts()
            return
        pts=np.asarray(self._all_points,dtype=float)
        self._tree=cKDTree(pts)
        self._indexed_keys=list(self._all_keys)
        self._indexed_count=len(self._all_keys)
        self.diagnostics.rebuilds+=1
        self._refresh_counts()

    def _refresh_counts(self):
        self.diagnostics.indexed_points=int(self._indexed_count)
        self.diagnostics.buffered_points=int(
            len(self._all_keys)-self._indexed_count
        )

    def nearest(self,q):
        q=np.asarray(q,dtype=float)
        if q.shape!=(self.dimension,):
            raise ValueError("query coordinate has incompatible shape.")
        if not self._all_keys:
            return None,None

        self.diagnostics.nearest_queries+=1
        best_key=None
        best_distance=np.inf

        if self._tree is not None and self._indexed_count:
            self.diagnostics.kd_queries+=1
            distance,index=self._tree.query(q,k=1)
            best_distance=float(distance)
            best_key=self._indexed_keys[int(index)]

        for key,point in zip(
            self._all_keys[self._indexed_count:],
            self._all_points[self._indexed_count:],
        ):
            self.diagnostics.buffer_distance_checks+=1
            distance=float(np.linalg.norm(q-point))
            if (
                distance<best_distance-1e-14
                or (
                    abs(distance-best_distance)<=1e-14
                    and (
                        best_key is None
                        or repr(key)<repr(best_key)
                    )
                )
            ):
                best_distance=distance
                best_key=key

        return best_key,float(best_distance)

    @property
    def size(self):
        return len(self._all_keys)


class IndexedTrackedMolecularDirectProviderV20(
    TrackedMolecularDirectProviderV19
):
    """v0.19 molecular provider with an exact buffered KD-tree anchor index."""

    def __init__(
        self,
        backend,
        geometry_map,
        *,
        overlap_engine=None,
        tracking=None,
        failure=None,
        rebuild_batch=32,
    ):
        kwargs={}
        if tracking is not None:
            kwargs["tracking"]=tracking
        if failure is not None:
            kwargs["failure"]=failure
        super().__init__(
            backend,
            geometry_map,
            overlap_engine=overlap_engine,
            **kwargs,
        )
        self.spatial_index=BufferedKDTreeIndexV20(
            geometry_map.nq,
            rebuild_batch=rebuild_batch,
        )

    def _after_cache_insert(self,key,q):
        self.spatial_index.insert(key,q)

    def _nearest_key(self,q):
        return self.spatial_index.nearest(q)

    def diagnostics_dict(self):
        out=super().diagnostics_dict()
        out["spatial_index"] = (
            self.spatial_index.diagnostics.as_dict()
        )
        return out
