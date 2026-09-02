import numpy as np

from gaussian_dynamics.dynamic_gauge_graph import AnalyticCI2DFrameProvider
from gaussian_dynamics.dynamic_graph_aims import DynamicGraphTBF
from gaussian_dynamics.born_huang_grid_v12 import build_born_huang_grid_2d
from gaussian_dynamics.error_controlled_sparse_dynamics_v17 import (
    ErrorControlledSparseSettingsV17,
    run_error_controlled_sparse_lvc_gaussians,
)


def _seed():
    return DynamicGraphTBF(
        uid=0,state=1,
        q=np.array([0.55,0.45]),
        p=np.array([0.6,0.8]),
        A=1.2*np.eye(2),
        node=("seed",0),
    )


def test_short_v17_run_audits_and_enriches():
    provider=AnalyticCI2DFrameProvider(
        nuclear_mass_au=20.0
    )
    grid=build_born_huang_grid_2d(
        grid_n=22,half_width=3.0,mass=20.0
    )

    out=run_error_controlled_sparse_lvc_gaussians(
        [_seed()],
        provider=provider,
        grid=grid,
        dt=2e-4,
        steps=8,
        settings=ErrorControlledSparseSettingsV17(
            defect_interval=2,
            enrich_relative_threshold=1e-5,
            prune_relative_threshold=1e-8,
            minimum_capture_fraction=1e-5,
            minimum_local_utility=0.0,
            min_basis=1,max_basis=3,
            minimum_adaptation_separation_steps=2,
            minimum_prune_age_steps=10,
            residual_shortlist=4,
            candidate_position_shifts=(0.0,0.05,-0.05),
            candidate_width_scales=(0.8,1.2),
            edge_enter_score=0.05,
            edge_exit_score=0.025,
            search_overlap_floor=1e-5,
            max_audit_S_error=0.05,
            max_audit_H_error=0.05,
            max_audit_Snuc_error=0.05,
            sparsity_audit_interval=2,
            condition_limit=1e8,
            hard_condition_limit=1e10,
        ),
        store_every=2,
    )

    assert out["audit_history"]
    assert out["complexity"]["dense_audits"]>=1

    events=[
        e for e in out["events"]
        if e["kind"]=="sparse_cost_aware_enrichment"
    ]
    assert events
    assert events[0]["relative_defect_after"]<events[0]["relative_defect_before"]

    norms=np.array([r["norm"] for r in out["records"]])
    assert np.max(np.abs(norms-1.0))<5e-4


def test_failed_audit_relaxes_graph_one_way():
    provider=AnalyticCI2DFrameProvider(
        nuclear_mass_au=20.0
    )
    grid=build_born_huang_grid_2d(
        grid_n=18,half_width=3.0,mass=20.0
    )

    basis=[
        _seed(),
        DynamicGraphTBF(
            uid=1,state=0,
            q=np.array([1.8,0.45]),
            p=np.array([-0.5,0.4]),
            A=1.2*np.eye(2),
            node=("seed",1),
        ),
    ]

    out=run_error_controlled_sparse_lvc_gaussians(
        basis,
        provider=provider,
        grid=grid,
        dt=1e-4,
        steps=1,
        settings=ErrorControlledSparseSettingsV17(
            defect_interval=100,
            enrich_relative_threshold=1.0,
            prune_relative_threshold=0.0,
            min_basis=2,max_basis=2,
            edge_enter_score=10.0,
            edge_exit_score=5.0,
            search_overlap_floor=1e-4,
            local_omitted_score_l2_budget=1e9,
            max_audit_S_error=1e-12,
            max_audit_H_error=1e-12,
            max_audit_Snuc_error=1e-12,
            audit_relaxation_factor=0.1,
            max_audit_relaxations=2,
            sparsity_audit_interval=1,
        ),
        store_every=1,
    )

    relax=[
        e for e in out["events"]
        if e["kind"]=="sparse_audit_relaxation"
    ]
    assert relax
    assert relax[-1]["new_enter_score"]<10.0
