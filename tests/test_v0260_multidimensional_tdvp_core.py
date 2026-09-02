from dataclasses import replace

import numpy as np
import pytest

from gaussian_dynamics.adaptive_multigaussian_tdvp_v252 import (
    ThawedGaussianSpinorStateV252,
    build_adaptive_gaussian_spinor_matrices_v252,
    build_adaptive_variational_metric_system_v252,
)
from gaussian_dynamics.multigaussian_tdvp_v251 import QuadraticSpinHamiltonianV251
from gaussian_dynamics.multidimensional_gaussian_tdvp_v260 import (
    DiagonalGaussianSpinorStateV260,
    MultidimensionalVariationalSettingsV260,
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
    return DiagonalGaussianSpinorStateV260(
        q=[[-0.25, 0.0]],
        p=[[3.0, 0.0]],
        widths=[[2.0, 2.0]],
        chirps=[[0.0, 0.0]],
        coefficients=[[1.0, 0.0]],
    ).normalized()


def test_v0260_one_dimensional_algebra_reduces_to_v0252():
    H0 = np.asarray([[0.01, 0.002j], [-0.002j, -0.01]])
    H1 = np.asarray([[0.003, 0.004], [0.004, -0.002]])
    H2 = np.asarray([[0.0005, 0.0], [0.0, 0.0007]])
    old_model = QuadraticSpinHamiltonianV251(900.0, H0, H1, H2).validate()
    new_model = QuadraticSpinHamiltonianNDV260(
        [[900.0]], H0, H1[None], H2[None, None]
    ).validate()
    old_state = ThawedGaussianSpinorStateV252(
        [-0.7, 0.8], [0.5, -0.3], [1.1, 0.8], [0.1, -0.05],
        [[0.7, 0.2j], [0.15 - 0.1j, 0.3]],
    ).normalized()
    new_state = DiagonalGaussianSpinorStateV260(
        old_state.q[:, None], old_state.p[:, None], old_state.widths[:, None],
        old_state.chirps[:, None], old_state.coefficients,
    ).validate(require_normalized=True)
    old_S, old_H = build_adaptive_gaussian_spinor_matrices_v252(old_state, old_model)
    new_S, new_H = build_multidimensional_gaussian_matrices_v260(new_state, new_model)
    old_metric = build_adaptive_variational_metric_system_v252(old_state, old_model)
    new_metric = build_multidimensional_metric_system_v260(new_state, new_model)
    assert np.max(np.abs(old_S - new_S)) < 3.0e-13
    assert np.max(np.abs(old_H - new_H)) < 3.0e-13
    assert np.max(np.abs(old_metric.metric - new_metric.metric)) < 3.0e-13
    assert np.max(np.abs(old_metric.rhs - new_metric.rhs)) < 3.0e-13
    assert np.max(np.abs(old_metric.velocity - new_metric.velocity)) < 3.0e-13


def test_v0260_metric_is_psd_and_svd_receipt_is_accurate():
    system = build_multidimensional_metric_system_v260(_state(), two_state_ci_soc_model_v260())
    assert np.min(np.linalg.eigvalsh(system.metric)) > -3.0e-10
    assert system.solve_receipt.linear_residual_relative < 3.0e-9
    assert system.solve_receipt.rank + system.solve_receipt.nullity == len(system.rhs)


def test_v0260_implicit_midpoint_is_signed_reversible():
    state = _state()
    model = two_state_ci_soc_model_v260()
    forward = multidimensional_implicit_midpoint_step_v260(state, model, 0.01)
    backward = multidimensional_implicit_midpoint_step_v260(forward.end, model, -0.01)
    assert np.max(
        np.abs(pack_multidimensional_parameters_v260(backward.end) - pack_multidimensional_parameters_v260(state))
    ) < 2.0e-8
    assert abs(forward.norm_change) < 1.0e-8
    assert abs(forward.energy_change_hartree) < 1.0e-8


def test_v0260_inactive_packet_shapes_are_frozen():
    state = DiagonalGaussianSpinorStateV260(
        q=[[-0.25, 0.0], [0.75, 0.0]],
        p=[[3.0, 0.0], [3.0, 0.0]],
        widths=[[2.0, 2.0], [2.0, 2.0]],
        chirps=[[0.0, 0.0], [0.0, 0.0]],
        coefficients=[[1.0, 0.0], [0.0, 0.0]],
    ).normalized()
    step = multidimensional_implicit_midpoint_step_v260(
        state, two_state_ci_soc_model_v260(), 0.01, active_shape_mask=[True, False]
    )
    assert np.array_equal(step.end.q[-1], step.start.q[-1])
    assert np.array_equal(step.end.p[-1], step.start.p[-1])
    assert np.array_equal(step.end.widths[-1], step.start.widths[-1])
    assert np.array_equal(step.end.chirps[-1], step.start.chirps[-1])


def test_v0260_full_correlated_width_request_fails_closed():
    with pytest.raises(ValueError, match="does not admit full correlated width matrices"):
        replace(
            MultidimensionalVariationalSettingsV260(),
            full_correlated_width_matrices=True,
        ).validate()


def test_v0260_anisotropic_width_rejects_general_rotation():
    state = DiagonalGaussianSpinorStateV260(
        [[0.0, 0.0]], [[0.0, 0.0]], [[1.0, 2.0]], [[0.0, 0.0]], [[1.0, 0.0]]
    ).normalized()
    root = 2.0**-0.5
    with pytest.raises(ValueError, match="not closed under general rotations"):
        state.coordinate_rotated([[root, -root], [root, root]])
