import numpy as np

from gaussian_dynamics.dynamic_gauge_graph import AnalyticCI2DFrameProvider
from gaussian_dynamics.dynamic_graph_aims import DynamicGraphTBF
from gaussian_dynamics.sparse_error_budget_v17 import (
    monotone_nonincreasing,
    score_threshold_snapshot_sweep,
    local_score_budget_snapshot_sweep,
    summarize_snapshot_convergence,
)


def _basis():
    return [
        DynamicGraphTBF(
            uid=i,
            state=i%2,
            q=np.array([0.2+0.8*i,0.15]),
            p=np.array([0.2*(-1)**i,0.1]),
            A=(1.0+0.1*i)*np.eye(2),
            node=("budget",i),
        )
        for i in range(4)
    ]


def test_monotone_helper():
    assert monotone_nonincreasing([3,2,2,1])
    assert not monotone_nonincreasing([3,2,2.1,1])


def test_snapshot_threshold_and_budget_sweeps_converge():
    provider=AnalyticCI2DFrameProvider(
        nuclear_mass_au=20.0
    )

    thresholds=score_threshold_snapshot_sweep(
        _basis(),provider,dt=0.005,
        enter_scores=(0.2,0.1,0.03,0.005),
        search_overlap_floor=1e-6,
    )
    budgets=local_score_budget_snapshot_sweep(
        _basis(),provider,dt=0.005,
        enter_score=0.2,
        budgets=(1e30,0.2,0.05,0.0),
        search_overlap_floor=1e-6,
    )
    summary=summarize_snapshot_convergence(
        thresholds,budgets
    )

    assert summary["threshold_S_monotone"]
    assert summary["threshold_H_monotone"]
    assert summary["budget_S_monotone"]
    assert summary["budget_H_monotone"]

    assert (
        thresholds[-1]["active_edges"]
        >=thresholds[0]["active_edges"]
    )
    assert (
        budgets[-1]["active_edges"]
        >=budgets[0]["active_edges"]
    )
