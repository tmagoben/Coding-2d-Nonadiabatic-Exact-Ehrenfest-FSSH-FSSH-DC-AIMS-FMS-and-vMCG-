import numpy as np
import pytest

from gaussian_dynamics.multidimensional_soc_v260 import (
    ExactGridSettingsV260,
    QuadraticSpinHamiltonianNDV260,
    UniformGrid2DV260,
    exact_grid_boundary_probability_v260,
    exact_grid_split_step_v260,
    initial_gaussian_spinor_2d_v260,
    kramers_doublet_ci_soc_model_v260,
    phase_aligned_grid_error_v260,
    run_exact_grid_ci_soc_v260,
    singlet_triplet_ci_soc_model_v260,
    two_state_ci_soc_model_v260,
)


def test_v0260_two_state_soc_gap_and_zero_soc_ci():
    active = two_state_ci_soc_model_v260(soc_scale=0.0025)
    zero = two_state_ci_soc_model_v260(soc_scale=0.0)
    assert np.diff(np.linalg.eigvalsh(active.H0))[0] == pytest.approx(0.005)
    assert np.diff(np.linalg.eigvalsh(zero.H0))[0] == pytest.approx(0.0)


def test_v0260_complete_doublet_has_kramers_pairs():
    model = kramers_doublet_ci_soc_model_v260()
    energies = np.linalg.eigvalsh(model.hamiltonian([0.3, -0.2]))
    assert model.nstate == 4
    assert energies[1] == pytest.approx(energies[0], abs=2.0e-14)
    assert energies[3] == pytest.approx(energies[2], abs=2.0e-14)


def test_v0260_complete_singlet_triplet_projectors_resolve_identity():
    model = singlet_triplet_ci_soc_model_v260()
    projectors = model.projectors
    assert model.nstate == 5
    assert np.trace(projectors["singlet"]).real == pytest.approx(2.0)
    assert np.trace(projectors["triplet"]).real == pytest.approx(3.0)
    assert np.allclose(projectors["singlet"] + projectors["triplet"], np.eye(5))


def test_v0260_model_derivative_matches_centered_difference():
    model = two_state_ci_soc_model_v260()
    q = np.asarray([0.31, -0.27])
    step = 1.0e-6
    observed = []
    for axis in range(2):
        delta = np.eye(2)[axis] * step
        observed.append((model.hamiltonian(q + delta) - model.hamiltonian(q - delta)) / (2 * step))
    assert np.max(np.abs(np.asarray(observed) - model.derivative(q))) < 2.0e-11


def test_v0260_coordinate_transform_preserves_potential():
    model = two_state_ci_soc_model_v260()
    angle = 0.4
    rotation = np.asarray([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])
    rotated = model.coordinate_rotated(rotation)
    q_new = np.asarray([0.2, -0.4])
    assert np.allclose(rotated.hamiltonian(q_new), model.hamiltonian(q_new @ rotation))


def test_v0260_invalid_model_contracts_fail_closed():
    model = two_state_ci_soc_model_v260()
    with pytest.raises(ValueError, match="positive definite"):
        QuadraticSpinHamiltonianNDV260(
            [[1.0, 0.0], [0.0, -1.0]], model.H0, model.H1, model.H2
        ).validate()
    with pytest.raises(ValueError, match="resolve the identity"):
        QuadraticSpinHamiltonianNDV260(
            model.mass_matrix_au,
            model.H0,
            model.H1,
            model.H2,
            projectors={"partial": np.diag([1.0, 0.0])},
        ).validate()


def test_v0260_exact_grid_is_unitary_reversible_and_boundary_clear():
    model = two_state_ci_soc_model_v260(mass_au=(50.0, 50.0))
    grid = UniformGrid2DV260.from_bounds((-6.0, 6.0), (-6.0, 6.0), (48, 48))
    psi0 = initial_gaussian_spinor_2d_v260(
        grid, [1.0, 0.0], center=(-0.25, 0.0), momentum=(3.0, 0.0), widths=(2.0, 2.0)
    )
    trajectory = run_exact_grid_ci_soc_v260(
        model, grid, psi0, settings=ExactGridSettingsV260(0.01, 6, 6)
    )
    reverse = trajectory.final_state
    for _ in range(6):
        reverse = exact_grid_split_step_v260(reverse, model, grid, -0.01)
    assert trajectory.maximum_norm_drift < 2.0e-12
    assert phase_aligned_grid_error_v260(psi0, reverse, grid) < 2.0e-11
    assert exact_grid_boundary_probability_v260(trajectory.final_state, grid) < 1.0e-12


def test_v0260_exact_grid_rejects_bad_grid_and_zero_timestep():
    with pytest.raises(ValueError, match="at least eight"):
        UniformGrid2DV260.from_bounds((-1, 1), (-1, 1), (7, 8))
    with pytest.raises(ValueError, match="nonzero"):
        ExactGridSettingsV260(dt_au=0.0).validate()
