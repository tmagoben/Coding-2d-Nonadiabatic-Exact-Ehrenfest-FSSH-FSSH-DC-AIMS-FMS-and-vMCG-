"""Physical scalar-gauge and packet-permutation controls for masked evolution."""
from dataclasses import replace

import numpy as np
import pytest

from gaussian_dynamics.complex_gauge_v21 import GaugeTransformedOperatorProviderV21, PhaseMixingGaugeV21
from gaussian_dynamics.correlated_gaussian_tdvp_v270 import CorrelatedGaussianSpinorStateV270
from gaussian_dynamics.curved_basis_adaptation_v293 import evaluate_curved_state_v293
from gaussian_dynamics.curved_basis_adaptation_v293 import coefficient_basis_v293
from gaussian_dynamics.molecular_active_v306 import molecular_active_midpoint_step_v306, activate_molecular_shapes_v306
from gaussian_dynamics.molecular_frame_bridge_v300 import MolecularFrameBridgeV300
from gaussian_dynamics.molecular_hamiltonian_v301 import build_molecular_grid_hamiltonian_v301
from gaussian_dynamics.molecular_spawn_v304 import _checkpoint_connection
from gaussian_dynamics.v302_benchmark import free_packet_fixture_v302
from gaussian_dynamics.v303_benchmark import soc_packet_fixture_v303


@pytest.fixture(scope="module")
def case():
    _, bridge, grid = free_packet_fixture_v302()
    state = CorrelatedGaussianSpinorStateV270(
        [[-.8], [.9]], [[.3], [-.2]], [[[1.4]], [[1.1]]],
        [[[.08]], [[-.03]]], [[.8], [.2j]],
    ).normalized()
    gauge = PhaseMixingGaugeV21(np.eye(1, dtype=complex), np.array([[.35]]), np.zeros(1))
    moving = MolecularFrameBridgeV300(GaugeTransformedOperatorProviderV21(bridge.provider, gauge), [0.0])
    moving_grid = build_molecular_grid_hamiltonian_v301(moving, grid.axes)
    transformed = replace(state, coefficients=np.array([gauge.matrix(q).conj().T @ c for q, c in zip(state.q, state.coefficients)]))
    return state, bridge, grid, gauge, transformed, moving, moving_grid


def test_masked_propagation_agrees_in_coordinate_dependent_gauge(case):
    state, bridge, grid, gauge, transformed, moving, moving_grid = case
    fixed_step = molecular_active_midpoint_step_v306(state, bridge, grid, .001, active_shape_mask=[True, False])
    moving_step = molecular_active_midpoint_step_v306(transformed, moving, moving_grid, .001, active_shape_mask=[True, False])
    a = evaluate_curved_state_v293(fixed_step.after, _checkpoint_connection(bridge), grid.quadrature)
    b = evaluate_curved_state_v293(moving_step.after, _checkpoint_connection(moving), moving_grid.quadrature)
    b = np.array([gauge.matrix(q) @ psi for q, psi in zip(grid.quadrature.points, b)])
    np.testing.assert_allclose(a, b, atol=5e-8, rtol=0)
    for result, original in ((fixed_step, state), (moving_step, transformed)):
        for name in ("q", "p", "width_matrices", "chirp_matrices"):
            np.testing.assert_array_equal(getattr(result.after, name)[-1], getattr(original, name)[-1])


def test_activation_decision_and_projected_shape_spectrum_are_gauge_covariant(case):
    state, bridge, grid, _, transformed, moving, moving_grid = case
    mask_a, a = activate_molecular_shapes_v306(state, bridge, grid, [True, False])
    mask_b, b = activate_molecular_shapes_v306(transformed, moving, moving_grid, [True, False])
    np.testing.assert_array_equal(mask_a, mask_b)
    assert mask_a.tolist() == [True, True]
    left, right = a["trials"][0], b["trials"][0]
    assert left["reason"] == right["reason"] == "activation_gates_passed"
    assert left["projected_shape_rank"] == right["projected_shape_rank"] == 4
    for key in ("coefficient_norm_squared", "projected_shape_condition", "shape_velocity_amplification"):
        assert left[key] == pytest.approx(right[key], rel=1e-4, abs=1e-9)


def test_packet_permutation_preserves_masked_physical_evolution(case):
    state, bridge, grid, *_ = case
    swapped = replace(state, **{name: getattr(state, name)[::-1] for name in ("q", "p", "width_matrices", "chirp_matrices", "coefficients")})
    a = molecular_active_midpoint_step_v306(state, bridge, grid, .001, active_shape_mask=[True, False])
    b = molecular_active_midpoint_step_v306(swapped, bridge, grid, .001, active_shape_mask=[False, True])
    for name in ("q", "p", "width_matrices", "chirp_matrices", "coefficients"):
        np.testing.assert_allclose(getattr(a.after, name), getattr(b.after, name)[::-1], atol=2e-8, rtol=0)


def test_complete_soc_spinor_frozen_basis_matches_independent_matrix_midpoint():
    _, bridge, grid, gauge, _ = soc_packet_fixture_v303()
    _, moving, moving_grid, _, _ = soc_packet_fixture_v303(moving=True)
    state = CorrelatedGaussianSpinorStateV270(
        [[-.8], [.9]], [[.3], [-.2]], [[[1.4]], [[1.1]]],
        [[[.08]], [[-.03]]], [[1., 0., 0., 0.], [0., .2j, 0., 0.]],
    ).normalized()
    phi = coefficient_basis_v293(state, _checkpoint_connection(bridge), grid.quadrature).reshape(-1, 8)
    weights = np.repeat(grid.quadrature.weights, state.nstate)
    gram = phi.conj().T @ (weights[:, None] * phi)
    operator = phi.conj().T @ (weights[:, None] * (grid.matrix @ phi))
    dt = .01
    expected = np.linalg.solve(gram+.5j*dt*operator, (gram-.5j*dt*operator) @ state.coefficients.ravel())
    fixed = molecular_active_midpoint_step_v306(state, bridge, grid, dt, active_shape_mask=[False, False])
    np.testing.assert_allclose(fixed.after.coefficients.ravel(), expected, atol=3e-10, rtol=0)
    assert len(fixed.receipt["active_parameter_indices"]) == 16
    transformed = replace(state, coefficients=np.array([gauge.matrix(q).conj().T @ c for q, c in zip(state.q, state.coefficients)]))
    rotated = molecular_active_midpoint_step_v306(transformed, moving, moving_grid, dt, active_shape_mask=[False, False])
    expected_rotated = np.array([gauge.matrix(q).conj().T @ c for q, c in zip(state.q, fixed.after.coefficients)])
    np.testing.assert_allclose(rotated.after.coefficients, expected_rotated, atol=2e-9, rtol=0)
