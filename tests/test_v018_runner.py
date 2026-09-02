import numpy as np

from gaussian_dynamics.dynamic_gauge_graph import AnalyticCI2DFrameProvider
from gaussian_dynamics.dynamic_graph_aims import DynamicGraphTBF
from gaussian_dynamics.born_huang_grid_v12 import build_born_huang_grid_2d
from gaussian_dynamics.convergence_complete_dynamics_v18 import (
    ConvergenceCompleteSettingsV18,
    run_convergence_complete_lvc_gaussians,
)


def _seed():
    return DynamicGraphTBF(
        uid=0,state=1,
        q=np.array([0.55,0.45]),
        p=np.array([0.6,0.8]),
        A=1.2*np.eye(2),
        node=("seed",0),
    )


def test_v18_uses_sampled_normal_audits_and_only_two_dense_sentinels():
    provider=AnalyticCI2DFrameProvider(
        nuclear_mass_au=20.0
    )
    grid=build_born_huang_grid_2d(
        grid_n=20,half_width=3.0,mass=20.0
    )

    out=run_convergence_complete_lvc_gaussians(
        [_seed()],
        provider=provider,
        grid=grid,
        dt=2e-4,
        steps=8,
        settings=ConvergenceCompleteSettingsV18(
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
            edge_enter_score=0.03,
            edge_exit_score=0.015,
            search_overlap_floor=1e-5,
            local_omitted_score_l2_budget=0.08,
            sampled_audit_interval=2,
            sampled_audit_priority_pairs=2,
            sampled_audit_random_pairs=2,
            sentinel_max_S_error=0.05,
            sentinel_max_H_error=0.05,
            sentinel_max_Snuc_error=0.05,
            candidate_batch_size=3,
            condition_limit=1e8,
            hard_condition_limit=1e10,
        ),
        store_every=2,
        return_snapshots=True,
    )

    c=out["complexity"]
    assert c["sentinel_dense_audits"]==2
    assert c["sampled_audits"]==4
    assert len(out["sentinel_audit_history"])==2
    assert len(out["sampled_audit_history"])==4
    assert len(out["snapshots"])==5
    assert c["candidate_batches"]>=1

    norms=np.asarray([r["norm"] for r in out["records"]])
    assert np.max(np.abs(norms-1.0))<5e-4


def test_v18_resolves_adaptive_cadence_in_physical_time():
    provider=AnalyticCI2DFrameProvider(
        nuclear_mass_au=20.0
    )
    grid=build_born_huang_grid_2d(
        grid_n=16,half_width=3.0,mass=20.0
    )

    settings=ConvergenceCompleteSettingsV18(
        defect_interval=99,
        defect_interval_time=0.002,
        enrich_relative_threshold=1.0,
        prune_relative_threshold=0.0,
        min_basis=1,max_basis=1,
        minimum_adaptation_separation_time=0.004,
        minimum_prune_age_time=0.006,
        sampled_audit_interval_time=0.004,
        cost_horizon_time=0.005,
        edge_enter_score=0.03,
        edge_exit_score=0.015,
        sentinel_max_S_error=0.1,
        sentinel_max_H_error=0.1,
        sentinel_max_Snuc_error=0.1,
    )

    out=run_convergence_complete_lvc_gaussians(
        [_seed()],
        provider=provider,
        grid=grid,
        dt=0.001,
        steps=4,
        settings=settings,
        store_every=4,
    )

    resolved=out["settings"]["resolved_control_steps"]
    assert resolved["defect_interval"]==2
    assert resolved["sampled_audit_interval"]==4
    assert resolved["adaptation_separation"]==4
    assert resolved["prune_age"]==6
    assert resolved["cost_horizon"]==5
