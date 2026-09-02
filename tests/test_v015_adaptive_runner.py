import numpy as np

from gaussian_dynamics.dynamic_gauge_graph import AnalyticCI2DFrameProvider
from gaussian_dynamics.dynamic_graph_aims import DynamicGraphTBF
from gaussian_dynamics.born_huang_grid_v12 import build_born_huang_grid_2d
from gaussian_dynamics.adaptive_defect_dynamics_v14 import (
    AdaptiveDefectSettings,
    run_time_adaptive_defect_lvc_gaussians,
)
from gaussian_dynamics.adaptive_defect_dynamics_v15 import (
    AdaptiveDefectSettingsV15,
    run_time_adaptive_cost_aware_lvc_gaussians,
    reduced_density_from_snuc,
)


def _basis():
    return [
        DynamicGraphTBF(
            uid=0,state=1,
            q=np.array([0.55,0.45]),
            p=np.array([0.6,0.8]),
            A=1.2*np.eye(2),
            node=("seed",0),
        ),
        DynamicGraphTBF(
            uid=1,state=1,
            q=np.array([0.35,0.45]),
            p=np.array([0.6,0.8]),
            A=1.8*np.eye(2),
            node=("seed",1),
        ),
    ]


def test_cached_v15_fixed_basis_matches_v14_propagation():
    provider=AnalyticCI2DFrameProvider(
        nuclear_mass_au=20.0
    )
    grid=build_born_huang_grid_2d(
        grid_n=20,half_width=3.0,mass=20.0
    )

    old=run_time_adaptive_defect_lvc_gaussians(
        _basis(),
        provider=provider,
        grid=grid,
        dt=2e-4,
        steps=6,
        settings=AdaptiveDefectSettings(
            defect_interval=100,
            enrich_relative_threshold=1.0,
            prune_relative_threshold=0.0,
            min_basis=2,max_basis=2,
            check_initial_defect=False,
        ),
        store_every=6,
    )

    new=run_time_adaptive_cost_aware_lvc_gaussians(
        _basis(),
        provider=provider,
        grid=grid,
        dt=2e-4,
        steps=6,
        settings=AdaptiveDefectSettingsV15(
            defect_interval=100,
            enrich_relative_threshold=1.0,
            prune_relative_threshold=0.0,
            min_basis=2,max_basis=2,
            check_initial_defect=False,
        ),
        store_every=6,
    )

    rho_old=old["records"][-1]["diabatic_populations"]
    rho_new=new["records"][-1]["diabatic_populations"]

    assert np.allclose(rho_new,rho_old,atol=2e-10)

    c=new["complexity"]
    assert c["factorization_reduction_fraction"]>0.7
    assert c["pair_factorizations"]>0
    assert c["v14_factorization_baseline"]>c["propagation_pair_factorizations"]


def test_short_v15_run_uses_cost_aware_incremental_enrichment():
    provider=AnalyticCI2DFrameProvider(
        nuclear_mass_au=20.0
    )
    grid=build_born_huang_grid_2d(
        grid_n=24,half_width=3.0,mass=20.0
    )

    out=run_time_adaptive_cost_aware_lvc_gaussians(
        [_basis()[0]],
        provider=provider,
        grid=grid,
        dt=2e-4,
        steps=8,
        settings=AdaptiveDefectSettingsV15(
            defect_interval=2,
            enrich_relative_threshold=1e-5,
            prune_relative_threshold=1e-8,
            minimum_capture_fraction=1e-5,
            minimum_cost_aware_utility=0.0,
            min_basis=1,max_basis=3,
            minimum_adaptation_separation_steps=2,
            minimum_prune_age_steps=10,
            residual_shortlist=4,
            candidate_position_shifts=(0.0,0.05,-0.05),
            candidate_width_scales=(0.8,1.2),
            condition_limit=1e8,
            hard_condition_limit=1e10,
            check_initial_defect=False,
        ),
        store_every=2,
    )

    events=[
        e for e in out["events"]
        if e["kind"]=="cost_aware_defect_enrichment"
    ]
    assert events
    assert events[0]["incremental_matrix_expansion"]
    assert events[0]["zero_coefficient_insertion"]
    assert events[0]["relative_defect_after"]<events[0]["relative_defect_before"]

    c=out["complexity"]
    assert c["incremental_expansions"]>=1
    assert c["cost_ranking_calls"]>=1
    assert c["candidate_pair_factorizations"]>0
    assert c["factorization_avoided"]>0

    norms=np.array([r["norm"] for r in out["records"]])
    assert np.max(np.abs(norms-1.0))<5e-5
