import numpy as np

from gaussian_dynamics.dynamic_gauge_graph import AnalyticCI2DFrameProvider
from gaussian_dynamics.dynamic_graph_aims import DynamicGraphTBF
from gaussian_dynamics.edge_importance_v17 import (
    EdgeImportanceSettingsV17,
    ErrorControlledGaussianLocalityGraphV17,
    safe_global_overlap_radius,
    pair_specific_overlap_upper_bound,
)
from gaussian_dynamics.pair_cache_v15 import GaussianPairCache


def _tbf(uid,q,p=(0.0,0.0),A=None):
    if A is None:
        A=np.eye(2)
    return DynamicGraphTBF(
        uid=uid,state=uid%2,
        q=np.asarray(q,float),
        p=np.asarray(p,float),
        A=np.asarray(A,float),
        node=("v17",uid),
    )


def test_safe_global_radius_rejects_pairs_below_overlap_floor():
    floor=1e-4
    a=0.8
    R=safe_global_overlap_radius(a,floor)

    bound=pair_specific_overlap_upper_bound(
        np.array([0.0,0.0]),
        np.array([R+0.5,0.0]),
        a,a,
    )
    assert bound<floor


def test_edge_score_can_retain_hamiltonian_relevant_pair_above_overlap_only_rule():
    provider=AnalyticCI2DFrameProvider(
        nuclear_mass_au=20.0
    )
    basis=[
        _tbf(0,[-1.0,0.0],[0.4,0.0]),
        _tbf(1,[1.0,0.0],[-0.4,0.0]),
    ]

    graph=ErrorControlledGaussianLocalityGraphV17(
        provider,
        dt=0.005,
        settings=EdgeImportanceSettingsV17(
            enter_score=0.05,
            exit_score=0.025,
            search_overlap_floor=1e-5,
            overlap_weight=0.0,
            hamiltonian_weight=1.0,
            time_connection_weight=0.0,
        ),
    )
    update=graph.update(basis)

    assert update.exact_pair_checks==1
    info=update.importance[(0,1)]
    assert info.hamiltonian_relative>0.0
    if info.score>=0.05:
        assert update.active_edges==((0,1),)


def test_score_hysteresis_retains_edge_until_exit_score():
    provider=AnalyticCI2DFrameProvider(
        nuclear_mass_au=20.0
    )
    settings=EdgeImportanceSettingsV17(
        enter_score=0.20,
        exit_score=0.08,
        search_overlap_floor=1e-6,
        overlap_weight=1.0,
        hamiltonian_weight=0.0,
        time_connection_weight=0.0,
    )
    graph=ErrorControlledGaussianLocalityGraphV17(
        provider,dt=0.005,settings=settings
    )

    # Unit-width zero-momentum overlap exp(-dq^2/4).
    u0=graph.update([
        _tbf(0,[0.1,0.1]),
        _tbf(1,[2.1,0.1]),
    ])
    assert len(u0.active_edges)==1

    # Between enter and exit.
    u1=graph.update([
        _tbf(0,[0.1,0.1]),
        _tbf(1,[2.9,0.1]),
    ])
    assert len(u1.active_edges)==1

    # Below exit.
    u2=graph.update([
        _tbf(0,[0.1,0.1]),
        _tbf(1,[3.5,0.1]),
    ])
    assert len(u2.active_edges)==0


def test_relaxation_is_one_sided_and_reproducible():
    provider=AnalyticCI2DFrameProvider(
        nuclear_mass_au=20.0
    )
    graph=ErrorControlledGaussianLocalityGraphV17(
        provider,dt=0.005,
        settings=EdgeImportanceSettingsV17(
            enter_score=0.04,
            exit_score=0.02,
            search_overlap_floor=1e-4,
        ),
    )
    graph.relax_scores(0.5)
    graph.relax_search_floor(0.1)

    assert np.isclose(graph.settings.enter_score,0.02)
    assert np.isclose(graph.settings.exit_score,0.01)
    assert np.isclose(graph.settings.search_overlap_floor,1e-5)


def test_global_local_importance_budget_promotes_edges():
    provider=AnalyticCI2DFrameProvider(
        nuclear_mass_au=20.0
    )
    basis=[
        _tbf(0,[0.2,0.1]),
        _tbf(1,[1.5,0.1]),
        _tbf(2,[2.8,0.1]),
    ]

    loose=ErrorControlledGaussianLocalityGraphV17(
        provider,dt=0.005,
        settings=EdgeImportanceSettingsV17(
            enter_score=10.0,
            exit_score=5.0,
            search_overlap_floor=1e-6,
            local_omitted_score_l2_budget=1e9,
        ),
    ).update(basis)

    budgeted=ErrorControlledGaussianLocalityGraphV17(
        provider,dt=0.005,
        settings=EdgeImportanceSettingsV17(
            enter_score=10.0,
            exit_score=5.0,
            search_overlap_floor=1e-6,
            local_omitted_score_l2_budget=0.05,
        ),
    ).update(basis)

    assert budgeted.budget_promoted_edges>0
    assert (
        budgeted.omitted_candidate_score_l2
        <=0.05+1e-12
    )
    assert len(budgeted.active_edges)>len(loose.active_edges)
