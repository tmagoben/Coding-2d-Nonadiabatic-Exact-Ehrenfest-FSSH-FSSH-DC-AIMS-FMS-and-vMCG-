import numpy as np

from gaussian_dynamics.benchmark_campaign import CIPassageConfig
from gaussian_dynamics.dynamic_gauge_graph import AnalyticCI2DFrameProvider
from gaussian_dynamics.born_huang_grid_v12 import build_born_huang_grid_2d
from gaussian_dynamics.exact_benchmark import localized_adiabatic_packet_2d
from gaussian_dynamics.initial_projection_v12 import (
    make_shifted_initial_gaussian_bank,
    project_grid_wavefunction_to_spinor_complete_basis,
)
from gaussian_dynamics.v12_benchmark import (
    V12AcceptanceThresholds,
    evaluate_v12_acceptance,
)


def test_shifted_bank_improves_initial_projection_over_one_center_gaussian():
    config=CIPassageConfig()
    provider=AnalyticCI2DFrameProvider(nuclear_mass_au=config.mass)
    grid=build_born_huang_grid_2d(
        grid_n=36,
        half_width=config.half_width,
        mass=config.mass,
    )
    psi=localized_adiabatic_packet_2d(
        grid.points,
        config.q_array(),
        config.p_array(),
        config.A_matrix(),
        state=config.state,
    )

    one=make_shifted_initial_gaussian_bank(
        config.q_array(),
        config.p_array(),
        config.A_matrix(),
        config.state,
        shifts=((0.0,0.0),),
    )
    five=make_shifted_initial_gaussian_bank(
        config.q_array(),
        config.p_array(),
        2.0*config.A_matrix(),
        config.state,
        shifts=(
            (0.0,0.0),
            (0.35,0.0),(-0.35,0.0),
            (0.0,0.35),(0.0,-0.35),
        ),
    )

    r1=project_grid_wavefunction_to_spinor_complete_basis(
        psi,grid.points,grid.dx,one,provider
    )
    r5=project_grid_wavefunction_to_spinor_complete_basis(
        psi,grid.points,grid.dx,five,provider
    )

    assert r5.fidelity > r1.fidelity
    assert r5.relative_residual < r1.relative_residual


def test_v12_acceptance_is_separately_sensitive_to_initial_and_dynamic_error():
    row={
        "initial_density_error":0.03,
        "projected_dynamics_density_error":2e-4,
        "target_density_error":0.04,
        "target_population_error":0.03,
        "coherence_phase_error":0.002,
        "max_norm_drift":1e-6,
        "max_condition_number":2e3,
    }

    ok=evaluate_v12_acceptance(
        row,V12AcceptanceThresholds()
    )
    assert ok["passed"]

    bad=dict(row)
    bad["projected_dynamics_density_error"]=0.01
    result=evaluate_v12_acceptance(
        bad,V12AcceptanceThresholds()
    )
    assert not result["passed"]
    assert result["checks"]["initial_density_representation"]
    assert not result["checks"]["projected_dynamics"]
