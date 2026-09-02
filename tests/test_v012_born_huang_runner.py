import numpy as np

from gaussian_dynamics.benchmark_campaign import CIPassageConfig
from gaussian_dynamics.dynamic_graph_aims import DynamicGraphTBF
from gaussian_dynamics.dynamic_gauge_graph import AnalyticCI2DFrameProvider
from gaussian_dynamics.born_huang_grid_v12 import (
    build_born_huang_grid_2d,
    born_huang_reduced_density,
)
from gaussian_dynamics.born_huang_dynamics_v12 import (
    run_born_huang_projected_gaussians,
)


def _seed(config):
    return DynamicGraphTBF(
        uid=0,
        state=config.state,
        q=config.q_array(),
        p=config.p_array(),
        A=config.A_matrix(),
        node=("seed",0),
    )


def test_born_huang_short_run_conserves_metric_norm():
    config=CIPassageConfig(
        q0=(0.55,0.45),
        p0=(0.6,0.8),
        A_diag=(1.2,1.2),
        mass=20.0,
        final_time=0.004,
        half_width=3.0,
    )
    provider=AnalyticCI2DFrameProvider(nuclear_mass_au=config.mass)
    grid=build_born_huang_grid_2d(
        grid_n=24,
        half_width=config.half_width,
        mass=config.mass,
    )

    out=run_born_huang_projected_gaussians(
        [_seed(config)],
        [1.0+0j],
        provider=provider,
        grid=grid,
        dt=0.001,
        steps=4,
        spawn_action_threshold=1e9,
        max_basis=1,
        store_every=1,
    )

    drift=max(abs(r["norm"]-1.0) for r in out["records"])
    assert drift < 2e-5

    rho=born_huang_reduced_density(
        out["final_coefficients"],
        out["final_basis"],
        grid,
        normalize=True,
    )
    assert np.isclose(np.trace(rho),1.0)
    assert np.min(np.linalg.eigvalsh(rho)) > -1e-10


def test_born_huang_short_spawning_run_remains_finite():
    config=CIPassageConfig(
        q0=(0.55,0.45),
        p0=(0.6,0.8),
        A_diag=(1.2,1.2),
        mass=20.0,
        final_time=0.004,
        half_width=3.0,
    )
    provider=AnalyticCI2DFrameProvider(nuclear_mass_au=config.mass)
    grid=build_born_huang_grid_2d(
        grid_n=20,
        half_width=config.half_width,
        mass=config.mass,
    )

    out=run_born_huang_projected_gaussians(
        [_seed(config)],
        [1.0+0j],
        provider=provider,
        grid=grid,
        dt=0.001,
        steps=4,
        spawn_action_threshold=1e-7,
        overlap_block=0.99999,
        child_overlap_block=0.99999,
        max_basis=3,
        children_per_event=1,
        position_shifts=(0.0,0.04,-0.04),
        width_scales=(0.7,1.0,1.5),
        store_every=1,
    )

    assert np.all(np.isfinite(out["final_coefficients"]))
    assert len(out["final_basis"]) >= 1
