from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class ElectronicCostEstimate:
    cost_units: float
    cache_hit: bool
    nearest_distance: float
    source: str

    def as_dict(self):
        return {
            "cost_units":float(self.cost_units),
            "cache_hit":bool(self.cache_hit),
            "nearest_distance":
                float(self.nearest_distance),
            "source":str(self.source),
        }


class UniformElectronicCostModel:
    """Constant provider-cost model, appropriate for the analytic LVC benchmark."""

    def __init__(self,cost_units=0.0,source="analytic_uniform"):
        self.cost_units=float(cost_units)
        self.source=str(source)

    def estimate(self,q):
        _=np.asarray(q,dtype=float)
        return ElectronicCostEstimate(
            cost_units=self.cost_units,
            cache_hit=True,
            nearest_distance=0.0,
            source=self.source,
        )

    def register(self,q):
        _=np.asarray(q,dtype=float)


class GeometryCacheElectronicCostModel:
    r"""Simple geometry-cache model for expensive electronic-structure providers.

    A candidate within `reuse_radius` of a registered geometry is charged
    `cached_cost_units`.  Otherwise it is charged `new_cost_units`.

    This does not attempt to predict actual CASSCF scaling from orbital counts.  It is
    a deliberately transparent bridge between the adaptive utility and measured
    electronic-structure cache behavior.
    """

    def __init__(
        self,
        cached_geometries=(),
        *,
        reuse_radius=0.05,
        cached_cost_units=0.05,
        new_cost_units=1.0,
        source="geometry_cache",
    ):
        self.reuse_radius=float(reuse_radius)
        self.cached_cost_units=float(cached_cost_units)
        self.new_cost_units=float(new_cost_units)
        self.source=str(source)

        if self.reuse_radius<0.0:
            raise ValueError("reuse_radius cannot be negative.")
        if self.cached_cost_units<0.0 or self.new_cost_units<0.0:
            raise ValueError("cost units cannot be negative.")

        self._geometries=[]
        for q in cached_geometries:
            self.register(q)

    @property
    def cached_geometries(self):
        return tuple(q.copy() for q in self._geometries)

    def register(self,q):
        q=np.asarray(q,dtype=float)
        if q.ndim!=1:
            raise ValueError("geometry coordinate must be a vector.")
        if self._geometries and q.shape!=self._geometries[0].shape:
            raise ValueError("registered geometry dimensions must match.")
        self._geometries.append(q.copy())

    def estimate(self,q):
        q=np.asarray(q,dtype=float)
        if q.ndim!=1:
            raise ValueError("geometry coordinate must be a vector.")

        if not self._geometries:
            return ElectronicCostEstimate(
                cost_units=self.new_cost_units,
                cache_hit=False,
                nearest_distance=float("inf"),
                source=self.source,
            )

        distances=np.array([
            np.linalg.norm(q-x)
            for x in self._geometries
        ],dtype=float)
        nearest=float(np.min(distances))
        hit=nearest<=self.reuse_radius

        return ElectronicCostEstimate(
            cost_units=(
                self.cached_cost_units
                if hit
                else self.new_cost_units
            ),
            cache_hit=bool(hit),
            nearest_distance=nearest,
            source=self.source,
        )
