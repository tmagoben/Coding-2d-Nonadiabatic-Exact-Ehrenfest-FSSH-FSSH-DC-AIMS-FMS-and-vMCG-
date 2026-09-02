import numpy as np
from dataclasses import dataclass

from gaussian_dynamics.dynamic_graph_aims import DynamicGraphTBF
from gaussian_dynamics.electronic_cost_v16 import (
    UniformElectronicCostModel,
    GeometryCacheElectronicCostModel,
)
from gaussian_dynamics.local_cost_aware_v16 import (
    estimate_local_sparse_incremental_cost,
    rank_local_sparse_candidates,
)


def _tbf(uid,q):
    return DynamicGraphTBF(
        uid=uid,state=0,
        q=np.asarray(q,float),
        p=np.zeros(2),
        A=np.eye(2),
        node=("t",uid),
    )


def test_geometry_cache_cost_distinguishes_hit_and_new_point():
    model=GeometryCacheElectronicCostModel(
        [[0.0,0.0]],
        reuse_radius=0.1,
        cached_cost_units=0.05,
        new_cost_units=2.0,
    )

    hit=model.estimate([0.05,0.0])
    miss=model.estimate([0.5,0.0])

    assert hit.cache_hit
    assert hit.cost_units==0.05
    assert not miss.cache_hit
    assert miss.cost_units==2.0


def test_local_cost_uses_degree_not_global_basis_size():
    basis=[
        _tbf(0,[0,0]),
        _tbf(1,[0.5,0]),
        _tbf(2,[5.0,0]),
        _tbf(3,[8.0,0]),
    ]
    candidate=_tbf(99,[0.2,0])

    cost=estimate_local_sparse_incremental_cost(
        candidate,
        basis,
        active_offdiagonal_edges=1,
        overlap_threshold=0.2,
        horizon_steps=5,
    )
    assert cost.local_degree<4
    assert cost.additional_pair_factorizations==2*5*(cost.local_degree+1)


@dataclass(frozen=True)
class _ResidualScore:
    candidate_index: int
    capture_fraction: float
    expanded_condition_number: float
    parent_uid: int
    target_state: int
    label: str


@dataclass(frozen=True)
class _Candidate:
    candidate: object


class _GaussianCandidateProxy:
    def __init__(self,tbf):
        self.tbf=tbf
    def to_tbf(self,uid,node_prefix="x"):
        return _tbf(uid,self.tbf.q)


def test_electronic_cache_miss_can_change_cost_aware_order():
    basis=[
        _tbf(0,[0,0]),
        _tbf(1,[0.5,0]),
    ]
    dynamic=[
        _Candidate(_GaussianCandidateProxy(_tbf(10,[0.05,0]))),
        _Candidate(_GaussianCandidateProxy(_tbf(11,[0.8,0]))),
    ]
    scores=[
        _ResidualScore(0,0.12,10.0,0,1,"cached"),
        _ResidualScore(1,0.13,10.0,0,1,"new"),
    ]
    model=GeometryCacheElectronicCostModel(
        [[0,0]],
        reuse_radius=0.1,
        cached_cost_units=0.0,
        new_cost_units=2.0,
    )

    ranked=rank_local_sparse_candidates(
        scores,
        dynamic,
        basis,
        active_offdiagonal_edges=1,
        overlap_threshold=0.05,
        current_condition=10.0,
        electronic_cost_model=model,
        electronic_cost_weight=1.0,
    )

    assert ranked[0].candidate_index==0
    assert ranked[0].electronic_cache_hit
