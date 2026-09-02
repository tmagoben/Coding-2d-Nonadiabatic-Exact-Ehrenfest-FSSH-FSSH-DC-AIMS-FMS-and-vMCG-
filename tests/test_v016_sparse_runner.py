import numpy as np

from gaussian_dynamics.dynamic_gauge_graph import AnalyticCI2DFrameProvider
from gaussian_dynamics.dynamic_graph_aims import DynamicGraphTBF
from gaussian_dynamics.born_huang_grid_v12 import build_born_huang_grid_2d
from gaussian_dynamics.sparse_adaptive_dynamics_v16 import (
    SparseAdaptiveSettingsV16,
    run_sparse_cost_aware_lvc_gaussians,
)


def _seed():
    return DynamicGraphTBF(
        uid=0,state=1,
        q=np.array([0.55,0.45]),
        p=np.array([0.6,0.8]),
        A=1.2*np.eye(2),
        node=("seed",0),
    )


def test_short_sparse_run_enriches_and_records_graph_sparsity():
    provider=AnalyticCI2DFrameProvider(
        nuclear_mass_au=20.0
    )
    grid=build_born_huang_grid_2d(
        grid_n=22,half_width=3.0,mass=20.0
    )

    out=run_sparse_cost_aware_lvc_gaussians(
        [_seed()],
        provider=provider,
        grid=grid,
        dt=2e-4,
        steps=8,
        settings=SparseAdaptiveSettingsV16(
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
            locality_enter_overlap=0.1,
            locality_exit_overlap=0.05,
            condition_limit=1e8,
            hard_condition_limit=1e10,
        ),
        store_every=2,
    )

    events=[
        e for e in out["events"]
        if e["kind"]=="sparse_cost_aware_enrichment"
    ]
    assert events
    assert events[0]["relative_defect_after"]<events[0]["relative_defect_before"]

    norms=np.array([r["norm"] for r in out["records"]])
    assert np.max(np.abs(norms-1.0))<2e-4

    c=out["complexity"]
    assert c["endpoint_graph_updates"]>0
    assert c["midpoint_graph_updates"]>0
    assert c["sparse_cayley_solves"]==8
