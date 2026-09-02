import numpy as np
import pytest

from gaussian_dynamics.benchmark_campaign import CIPassageConfig
from gaussian_dynamics.dynamic_gauge_graph import AnalyticCI2DFrameProvider
from gaussian_dynamics.dynamic_graph_aims import DynamicGraphTBF
from gaussian_dynamics.born_huang_grid_v12 import build_born_huang_grid_2d
from gaussian_dynamics.adaptive_defect_dynamics_v14 import (
    AdaptiveDefectSettings,
    run_time_adaptive_defect_lvc_gaussians,
)


def _seed():
    return DynamicGraphTBF(
        uid=0,state=1,
        q=np.array([0.55,0.45]),
        p=np.array([0.6,0.8]),
        A=1.2*np.eye(2),
        node=("seed",0),
    )


def test_adaptive_settings_require_hysteresis():
    with pytest.raises(ValueError):
        AdaptiveDefectSettings(
            enrich_relative_threshold=0.01,
            prune_relative_threshold=0.02,
        ).validate()


def test_short_adaptive_run_enriches_and_records_complexity():
    provider=AnalyticCI2DFrameProvider(
        nuclear_mass_au=20.0
    )
    grid=build_born_huang_grid_2d(
        grid_n=24,
        half_width=3.0,
        mass=20.0,
    )

    settings=AdaptiveDefectSettings(
        defect_interval=2,
        enrich_relative_threshold=1e-5,
        prune_relative_threshold=1e-8,
        minimum_capture_fraction=1e-5,
        min_basis=1,
        max_basis=3,
        minimum_adaptation_separation_steps=2,
        minimum_prune_age_steps=10,
        candidate_position_shifts=(0.0,0.05,-0.05),
        candidate_width_scales=(0.8,1.2),
        condition_limit=1e8,
        hard_condition_limit=1e10,
    )

    out=run_time_adaptive_defect_lvc_gaussians(
        [_seed()],
        provider=provider,
        grid=grid,
        dt=2e-4,
        steps=8,
        settings=settings,
        store_every=2,
    )

    assert any(
        e["kind"]=="defect_enrichment"
        for e in out["events"]
    )
    assert len(out["final_basis"])>=2

    norms=np.array([r["norm"] for r in out["records"]])
    assert np.max(np.abs(norms-1.0))<5e-5

    complexity=out["complexity"]
    assert complexity["defect_evaluations"]>=1
    assert complexity["candidate_ranking_calls"]>=1
    assert complexity["candidate_count_scored"]>0
    assert complexity["peak_basis_size"]>=2
    assert complexity["cayley_solve_calls"]==8
