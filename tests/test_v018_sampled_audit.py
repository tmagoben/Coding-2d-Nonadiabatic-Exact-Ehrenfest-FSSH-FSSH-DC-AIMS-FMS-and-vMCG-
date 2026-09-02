import numpy as np

from gaussian_dynamics.dynamic_gauge_graph import AnalyticCI2DFrameProvider
from gaussian_dynamics.dynamic_graph_aims import DynamicGraphTBF
from gaussian_dynamics.edge_importance_v17 import (
    EdgeImportanceSettingsV17,
    ErrorControlledGaussianLocalityGraphV17,
)
from gaussian_dynamics.sampled_sparse_audit_v18 import (
    sampled_omitted_edge_audit_v18,
)


def _basis(n=8):
    out=[]
    for i in range(n):
        out.append(
            DynamicGraphTBF(
                uid=i,
                state=i%2,
                q=np.array([0.35+0.75*i,0.2]),
                p=np.array([0.05*(-1)**i,0.02]),
                A=np.eye(2),
                node=("audit",i),
            )
        )
    return out


def test_sampled_audit_is_deterministic():
    provider=AnalyticCI2DFrameProvider(
        nuclear_mass_au=20.0
    )
    settings=EdgeImportanceSettingsV17(
        enter_score=0.3,
        exit_score=0.15,
        search_overlap_floor=1e-4,
        local_omitted_score_l2_budget=1e9,
    )
    graph=ErrorControlledGaussianLocalityGraphV17(
        provider,0.005,settings
    )
    basis=_basis()
    update=graph.update(basis)

    a=sampled_omitted_edge_audit_v18(
        basis,provider,0.005,update,graph.settings,
        step=20,priority_count=3,random_count=3,
        seed=1234,
    )
    b=sampled_omitted_edge_audit_v18(
        basis,provider,0.005,update,graph.settings,
        step=20,priority_count=3,random_count=3,
        seed=1234,
    )

    assert a.sampled_pairs==b.sampled_pairs
    assert np.isclose(a.maximum_score,b.maximum_score)
    assert len(a.sampled_pairs)<=6


def test_sampled_audit_detects_missed_high_score_edge():
    provider=AnalyticCI2DFrameProvider(
        nuclear_mass_au=20.0
    )
    basis=[
        DynamicGraphTBF(
            uid=0,state=0,
            q=np.array([0.2,0.2]),
            p=np.array([0.6,0.0]),
            A=np.eye(2),
            node=("a",0),
        ),
        DynamicGraphTBF(
            uid=1,state=1,
            q=np.array([1.0,0.2]),
            p=np.array([-0.6,0.0]),
            A=np.eye(2),
            node=("a",1),
        ),
    ]

    # Force omission with an unrealistically high graph threshold.
    settings=EdgeImportanceSettingsV17(
        enter_score=10.0,
        exit_score=5.0,
        search_overlap_floor=1e-6,
        local_omitted_score_l2_budget=1e9,
    )
    graph=ErrorControlledGaussianLocalityGraphV17(
        provider,0.005,settings
    )
    update=graph.update(basis)
    assert update.active_edges==()

    # Audit against a much stricter violation threshold than the graph's artificial
    # enter score: any physically nonzero pair should be flagged.
    a=sampled_omitted_edge_audit_v18(
        basis,provider,0.005,update,graph.settings,
        step=0,priority_count=1,random_count=0,
        seed=7,violation_factor=1e-6,
    )
    assert not a.passed
    assert a.violation_count==1
