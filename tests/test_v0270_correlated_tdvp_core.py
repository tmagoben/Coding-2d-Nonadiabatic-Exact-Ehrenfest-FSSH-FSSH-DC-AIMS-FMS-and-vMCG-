from dataclasses import replace

import numpy as np
import pytest

from gaussian_dynamics.correlated_gaussian_tdvp_v270 import (
    CorrelatedGaussianSpinorStateV270,
    build_correlated_gaussian_matrices_v270,
    build_correlated_metric_system_v270,
    correlated_implicit_midpoint_step_v270,
    exp_frechet_symmetric_v270,
    exp_symmetric_v270,
    log_spd_v270,
    pack_correlated_parameters_v270,
    rotate_correlated_velocity_v270,
    smat_v270,
    state_from_correlated_parameters_v270,
    symmetric_basis_v270,
    svec_v270,
)
from gaussian_dynamics.multidimensional_gaussian_tdvp_v260 import (
    DiagonalGaussianSpinorStateV260,
    build_multidimensional_gaussian_matrices_v260,
    build_multidimensional_metric_system_v260,
    multidimensional_implicit_midpoint_step_v260,
    pack_multidimensional_parameters_v260,
)
from gaussian_dynamics.multidimensional_soc_v260 import (
    QuadraticSpinHamiltonianNDV260,
    two_state_ci_soc_model_v260,
)


def _state():
    return CorrelatedGaussianSpinorStateV270(
        [[-0.25, 0.1]], [[3.0, 0.2]],
        [[[1.7, 0.25], [0.25, 2.6]]],
        [[[0.02, 0.03], [0.03, -0.01]]],
        [[1.0, 0.0]],
    ).normalized()


def test_v0270_svec_and_log_width_coordinates_roundtrip():
    matrix = np.asarray([[0.3, -0.2], [-0.2, 0.7]])
    assert np.max(np.abs(smat_v270(svec_v270(matrix), 2) - matrix)) < 3.0e-14
    width = exp_symmetric_v270(matrix)
    assert np.max(np.abs(log_spd_v270(width) - matrix)) < 3.0e-14
    basis = symmetric_basis_v270(2)
    gram = np.asarray([[np.vdot(a, b) for b in basis] for a in basis])
    assert np.max(np.abs(gram - np.eye(3))) < 3.0e-14


def test_v0270_exp_frechet_matches_centered_difference():
    matrix = np.asarray([[0.31, -0.18], [-0.18, -0.27]])
    direction = np.asarray([[0.23, -0.11], [-0.11, 0.36]])
    epsilon = 2.0e-6
    finite = (
        exp_symmetric_v270(matrix + epsilon * direction)
        - exp_symmetric_v270(matrix - epsilon * direction)
    ) / (2.0 * epsilon)
    analytic = exp_frechet_symmetric_v270(matrix, direction)
    assert np.max(np.abs(finite - analytic)) < 2.0e-8


def test_v0270_pack_roundtrip_preserves_full_matrices():
    state = _state()
    recovered = state_from_correlated_parameters_v270(
        pack_correlated_parameters_v270(state), ngaussian=state.ngaussian,
        ndim=state.ndim, nstate=state.nstate, time_au=state.time_au,
    )
    assert np.max(
        np.abs(pack_correlated_parameters_v270(recovered) - pack_correlated_parameters_v270(state))
    ) < 3.0e-13
    assert recovered.width_matrices[0, 0, 1] != 0.0
    assert recovered.chirp_matrices[0, 0, 1] != 0.0


def test_v0270_analytic_matrices_are_hermitian():
    overlap, hamiltonian = build_correlated_gaussian_matrices_v270(
        _state(), two_state_ci_soc_model_v260()
    )
    assert np.max(np.abs(overlap - overlap.conj().T)) < 3.0e-12
    assert np.max(np.abs(hamiltonian - hamiltonian.conj().T)) < 3.0e-12


def test_v0270_one_dimensional_dynamics_reduce_exactly_to_v0260():
    H0 = np.asarray([[0.01, 0.002j], [-0.002j, -0.01]])
    H1 = np.asarray([[[0.003, 0.004], [0.004, -0.002]]])
    H2 = np.asarray([[[[0.0005, 0.0], [0.0, 0.0007]]]])
    model = QuadraticSpinHamiltonianNDV260([[900.0]], H0, H1, H2).validate()
    old = DiagonalGaussianSpinorStateV260(
        [[-0.7], [0.8]], [[0.5], [-0.3]], [[1.1], [0.8]], [[0.1], [-0.05]],
        [[0.7, 0.2j], [0.15 - 0.1j, 0.3]],
    ).normalized()
    new = CorrelatedGaussianSpinorStateV270.from_diagonal_v260(old)
    old_matrices = build_multidimensional_gaussian_matrices_v260(old, model)
    new_matrices = build_correlated_gaussian_matrices_v270(new, model)
    old_metric = build_multidimensional_metric_system_v260(old, model)
    new_metric = build_correlated_metric_system_v270(new, model)
    assert np.max(np.abs(old_matrices[0] - new_matrices[0])) < 3.0e-12
    assert np.max(np.abs(old_matrices[1] - new_matrices[1])) < 3.0e-12
    assert np.max(np.abs(old_metric.metric - new_metric.metric)) < 3.0e-12
    assert np.max(np.abs(old_metric.velocity - new_metric.velocity)) < 3.0e-12
    old_step = multidimensional_implicit_midpoint_step_v260(old, model, 0.001)
    new_step = correlated_implicit_midpoint_step_v270(new, model, 0.001)
    assert np.max(
        np.abs(
            pack_multidimensional_parameters_v260(old_step.end)
            - pack_correlated_parameters_v270(new_step.end)
        )
    ) < 3.0e-11


def test_v0270_metric_and_midpoint_are_arbitrary_rotation_covariant():
    state = _state()
    model = two_state_ci_soc_model_v260()
    angle = 0.371
    rotation = np.asarray(
        [[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]]
    )
    base = build_correlated_metric_system_v270(state, model)
    transformed_state = state.coordinate_rotated(rotation)
    transformed = build_correlated_metric_system_v270(
        transformed_state, model.coordinate_rotated(rotation)
    )
    expected_velocity = rotate_correlated_velocity_v270(state, base.velocity, rotation)
    assert np.max(np.abs(expected_velocity - transformed.velocity)) < 3.0e-10
    base_step = correlated_implicit_midpoint_step_v270(state, model, 0.002)
    transformed_step = correlated_implicit_midpoint_step_v270(
        transformed_state, model.coordinate_rotated(rotation), 0.002
    )
    assert np.max(
        np.abs(
            pack_correlated_parameters_v270(base_step.end.coordinate_rotated(rotation))
            - pack_correlated_parameters_v270(transformed_step.end)
        )
    ) < 3.0e-8


def test_v0270_midpoint_is_signed_reversible():
    state = _state()
    model = two_state_ci_soc_model_v260()
    forward = correlated_implicit_midpoint_step_v270(state, model, 0.002)
    backward = correlated_implicit_midpoint_step_v270(forward.end, model, -0.002)
    assert np.max(
        np.abs(pack_correlated_parameters_v270(backward.end) - pack_correlated_parameters_v270(state))
    ) < 3.0e-8
    assert abs(forward.norm_change) < 1.0e-8
    assert abs(forward.energy_change_hartree) < 1.0e-8


def test_v0270_inactive_full_matrix_shape_is_bitwise_frozen():
    state = CorrelatedGaussianSpinorStateV270(
        [[-0.25, 0.1], [0.8, -0.2]], [[3.0, 0.2], [2.0, -0.1]],
        [[[1.7, 0.25], [0.25, 2.6]], [[1.2, -0.2], [-0.2, 2.1]]],
        [[[0.02, 0.03], [0.03, -0.01]], [[-0.03, 0.02], [0.02, 0.04]]],
        [[1.0, 0.0], [0.0, 0.0]],
    ).normalized()
    step = correlated_implicit_midpoint_step_v270(
        state, two_state_ci_soc_model_v260(), 0.001, active_shape_mask=[True, False]
    )
    assert np.array_equal(step.start.q[-1], step.end.q[-1])
    assert np.array_equal(step.start.p[-1], step.end.p[-1])
    assert np.array_equal(step.start.width_matrices[-1], step.end.width_matrices[-1])
    assert np.array_equal(step.start.chirp_matrices[-1], step.end.chirp_matrices[-1])


def test_v0270_invalid_width_and_coordinate_transform_fail_closed():
    state = _state()
    with pytest.raises(ValueError, match="positive definite"):
        replace(state, width_matrices=np.asarray([[[1.0, 0.0], [0.0, -1.0]]])).validate()
    with pytest.raises(ValueError, match="orthogonal"):
        state.coordinate_rotated([[1.0, 0.2], [0.0, 1.0]])
    with pytest.raises(ValueError, match="cannot be represented"):
        state.to_diagonal_v260()
