import numpy as np

from gaussian_dynamics.benchmark_campaign import CIPassageConfig
from gaussian_dynamics.dynamic_gauge_graph import AnalyticCI2DFrameProvider
from gaussian_dynamics.dynamic_graph_aims import DynamicGraphTBF
from gaussian_dynamics.born_huang_grid_v12 import build_born_huang_grid_2d
from gaussian_dynamics.exact_benchmark import localized_adiabatic_packet_2d
from gaussian_dynamics.residual_basis_v13 import (
    cartesian_offsets_2d,
    generate_gaussian_dictionary,
    build_residual_greedy_basis,
    prepare_gaussian_dictionary,
    build_residual_greedy_basis_prepared,
)


def _problem(grid_n=32):
    c=CIPassageConfig()
    provider=AnalyticCI2DFrameProvider(nuclear_mass_au=c.mass)
    grid=build_born_huang_grid_2d(
        grid_n=grid_n,
        half_width=c.half_width,
        mass=c.mass,
    )
    psi=localized_adiabatic_packet_2d(
        grid.points,
        c.q_array(),
        c.p_array(),
        c.A_matrix(),
        state=c.state,
    )
    seed=DynamicGraphTBF(
        uid=0,
        state=c.state,
        q=c.q_array(),
        p=c.p_array(),
        A=c.A_matrix(),
        node=("seed",0),
    )
    candidates=generate_gaussian_dictionary(
        c.q_array(),
        c.p_array(),
        c.A_matrix(),
        c.state,
        cartesian_offsets_2d(radius=0.8,spacing=0.4),
        width_scales=(1.0,2.0,4.0),
    )
    return c,provider,grid,psi,seed,candidates


def test_residual_greedy_projection_improves_monotonically():
    _,provider,grid,psi,seed,candidates=_problem()

    result=build_residual_greedy_basis(
        psi,
        grid.points,
        grid.dx,
        provider,
        [seed],
        candidates,
        max_basis=5,
        top_k_density_screen=1,
        density_screen=False,
        condition_limit=1e6,
    )

    residuals=[step.relative_residual for step in result.history]
    assert len(residuals)==4
    assert all(
        residuals[i+1] < residuals[i]
        for i in range(len(residuals)-1)
    )
    assert result.final_relative_residual < 0.5


def test_predicted_gain_matches_actual_one_step_projection_gain():
    _,provider,grid,psi,seed,candidates=_problem()

    result=build_residual_greedy_basis(
        psi,
        grid.points,
        grid.dx,
        provider,
        [seed],
        candidates,
        max_basis=2,
        top_k_density_screen=1,
        density_screen=False,
        condition_limit=1e6,
    )

    step=result.history[0]
    target_norm=float(np.sum(np.abs(psi)**2)*grid.area)

    assert np.isclose(
        step.actual_residual_reduction,
        step.predicted_gain/target_norm,
        rtol=5e-5,
        atol=5e-7,
    )


def test_density_screened_residual_builder_is_deterministic():
    _,provider,grid,psi,seed,candidates=_problem()

    kwargs=dict(
        target_psi=psi,
        points=grid.points,
        dx=grid.dx,
        provider=provider,
        seed_basis=[seed],
        candidates=candidates,
        max_basis=4,
        top_k_density_screen=8,
        density_screen=True,
        condition_limit=1e6,
    )

    a=build_residual_greedy_basis(**kwargs)
    b=build_residual_greedy_basis(**kwargs)

    assert a.candidate_indices==b.candidate_indices
    assert np.isclose(
        a.final_relative_residual,
        b.final_relative_residual,
    )
    assert np.allclose(
        a.projection.coefficients,
        b.projection.coefficients,
    )


def test_prepared_residual_builder_matches_slow_greedy_choice():
    _,provider,grid,psi,seed,candidates=_problem(grid_n=28)

    slow=build_residual_greedy_basis(
        psi,grid.points,grid.dx,provider,[seed],candidates,
        max_basis=3,
        top_k_density_screen=1,
        density_screen=False,
        condition_limit=1e6,
    )

    prepared=prepare_gaussian_dictionary(
        candidates,grid.points,grid.dx
    )
    fast=build_residual_greedy_basis_prepared(
        psi,grid.points,grid.dx,provider,[seed],prepared,
        max_basis=3,
        top_k_density_screen=1,
        density_screen=False,
        condition_limit=1e6,
    )

    assert slow.candidate_indices==fast.candidate_indices
    assert np.isclose(
        slow.final_relative_residual,
        fast.final_relative_residual,
        rtol=2e-5,
        atol=2e-6,
    )
