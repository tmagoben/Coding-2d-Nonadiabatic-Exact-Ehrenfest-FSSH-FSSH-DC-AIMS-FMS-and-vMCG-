from dataclasses import replace
from types import SimpleNamespace

import numpy as np
import pytest

from gaussian_dynamics.analytic_soc_models_v220 import (
    AnalyticDoubletSOCProviderV220,
    AnalyticSingletTripletSOCProviderV220,
    DoubletSOCConfigV220,
    SingletTripletSOCConfigV220,
)
from gaussian_dynamics.complex_gauge_v21 import (
    GaugeTransformedOperatorProviderV21,
    PhaseMixingGaugeV21,
    random_unitary_v21,
)
from gaussian_dynamics.electronic_operator_v21 import (
    ElectronicOperatorPointV21,
    ElectronicOperatorSnapshotV21,
)
from gaussian_dynamics.variational_soc_dynamics_v250 import (
    ELECTRONIC_INTEGRATOR_V250,
    GENERAL_TDVP_INTEGRATOR_V250,
    POLAR_ALGORITHM_V250,
    CanonicalVariationalSOCStateV250,
    VariationalSOCIntegratorSettingsV250,
    reverse_variational_soc_trajectory_v250,
    run_symmetric_variational_soc_dynamics_v250,
    symmetric_variational_soc_step_v250,
)


def _coefficients():
    result = np.asarray(
        [0.67 + 0.11j, 0.19 - 0.28j, 0.41 + 0.17j, -0.09j],
        dtype=complex,
    )
    return result / np.linalg.norm(result)


def _state(q=-0.2, p=1.4):
    return CanonicalVariationalSOCStateV250(
        np.asarray([q]), np.asarray([p]), _coefficients()
    ).validate()


def _gauge():
    return PhaseMixingGaugeV21(
        random_unitary_v21(4, 25001),
        np.asarray([[0.17], [-0.11], [0.23], [-0.07]]),
        np.asarray([0.20, -0.30, 0.10, 0.40]),
    )


def test_scope_freezes_restricted_tdvp_and_recommends_midpoint_for_full_tdvp():
    settings = VariationalSOCIntegratorSettingsV250().validate()
    assert settings.electronic_integrator == ELECTRONIC_INTEGRATOR_V250
    assert settings.polar_algorithm == POLAR_ALGORITHM_V250
    assert "implicit midpoint" in GENERAL_TDVP_INTEGRATOR_V250

    with pytest.raises(ValueError, match="does not admit full multi-Gaussian TDVP"):
        replace(settings, full_multi_gaussian_tdvp=True).validate()
    with pytest.raises(ValueError, match="adaptive Gaussian-width"):
        replace(settings, adaptive_gaussian_widths=True).validate()
    with pytest.raises(ValueError, match="coordinate-dependent generalized mass"):
        replace(settings, coordinate_dependent_mass=True).validate()
    with pytest.raises(ValueError, match="SVD-polar algorithm is frozen"):
        replace(settings, polar_algorithm="raw overlap").validate()


def test_state_normalization_is_explicit_and_zero_vector_is_rejected():
    raw = CanonicalVariationalSOCStateV250(
        np.asarray([0.1]), np.asarray([0.2]), 2.0 * _coefficients()
    )
    with pytest.raises(ValueError, match="unit norm"):
        raw.validate()
    assert raw.normalized().electronic_norm == pytest.approx(1.0, abs=2.0e-16)
    with pytest.raises(ValueError, match="cannot normalize"):
        CanonicalVariationalSOCStateV250(
            np.asarray([0.1]), np.asarray([0.2]), np.zeros(4, dtype=complex)
        ).normalized()


def test_one_step_is_bound_to_verlet_strang_and_svd_polar_data():
    receipt = symmetric_variational_soc_step_v250(
        _state(), AnalyticDoubletSOCProviderV220(), 0.2
    )
    assert receipt.end.electronic_norm == pytest.approx(1.0, abs=2.0e-15)
    assert np.allclose(receipt.singular_values, 1.0)
    assert np.allclose(receipt.transport_end_to_start, np.eye(4))
    assert receipt.transport_metrics["trajectory_ready"] is True
    assert receipt.transport_policy["minimum_retained_singular_value"] == 0.9
    assert np.array_equal(receipt.mass_matrix_au, receipt.mass_matrix_end_au)
    assert receipt.as_dict()["transport_policy"] == receipt.transport_policy


@pytest.mark.parametrize(
    "provider",
    [AnalyticSingletTripletSOCProviderV220(), AnalyticDoubletSOCProviderV220()],
)
def test_forward_then_signed_adjoint_returns_initial_state(provider):
    initial = _state(q=-0.31, p=2.1)
    forward = run_symmetric_variational_soc_dynamics_v250(
        initial, provider, dt_au=0.3, steps=20
    )
    reverse = reverse_variational_soc_trajectory_v250(forward, provider)

    assert np.allclose(reverse.final_state.q, initial.q, atol=2.0e-14)
    assert np.allclose(reverse.final_state.p, initial.p, atol=2.0e-14)
    assert np.allclose(
        reverse.final_state.electronic_coefficients,
        initial.electronic_coefficients,
        atol=8.0e-14,
    )
    assert reverse.final_state.time_au == pytest.approx(0.0, abs=2.0e-14)


def test_coordinate_dependent_complex_gauge_covariance_is_complete():
    initial = _state(q=-0.27, p=3.2)
    base_provider = AnalyticDoubletSOCProviderV220()
    gauge = _gauge()
    gauge_provider = GaugeTransformedOperatorProviderV21(base_provider, gauge)
    transformed_initial = CanonicalVariationalSOCStateV250(
        initial.q,
        initial.p,
        gauge.matrix(initial.q).conj().T @ initial.electronic_coefficients,
    ).validate()

    base = run_symmetric_variational_soc_dynamics_v250(
        initial, base_provider, dt_au=0.4, steps=25
    )
    transformed = run_symmetric_variational_soc_dynamics_v250(
        transformed_initial, gauge_provider, dt_au=0.4, steps=25
    )
    expected = (
        gauge.matrix(transformed.final_state.q).conj().T
        @ base.final_state.electronic_coefficients
    )

    assert np.allclose(transformed.final_state.q, base.final_state.q, atol=2.0e-13)
    assert np.allclose(transformed.final_state.p, base.final_state.p, atol=2.0e-13)
    assert np.allclose(
        transformed.final_state.electronic_coefficients, expected, atol=2.0e-13
    )
    assert transformed.maximum_norm_drift < 8.0e-15


class _ScaledOverlapProvider:
    def __init__(self, singular_value):
        self.base = AnalyticDoubletSOCProviderV220()
        self.singular_value = float(singular_value)

    def evaluate_snapshot(self, q):
        return self.base.evaluate_snapshot(q)

    def snapshot_overlap(self, left, right):
        return self.singular_value * np.eye(4, dtype=complex)


def test_raw_contraction_is_retained_but_only_its_polar_factor_propagates():
    trajectory = run_symmetric_variational_soc_dynamics_v250(
        _state(), _ScaledOverlapProvider(0.97), dt_au=0.2, steps=5
    )
    for receipt in trajectory.steps:
        assert np.allclose(receipt.singular_values, 0.97)
        assert np.allclose(receipt.overlap_start_end, 0.97 * np.eye(4))
        assert np.allclose(receipt.transport_end_to_start, np.eye(4))
    assert trajectory.maximum_norm_drift < 3.0e-15


@pytest.mark.parametrize(
    "singular_value, match",
    [(1.01, "physically inconsistent"), (0.5, "not trajectory ready")],
)
def test_spectral_expansion_and_manifold_loss_fail_closed(singular_value, match):
    with pytest.raises(ValueError, match=match):
        symmetric_variational_soc_step_v250(
            _state(), _ScaledOverlapProvider(singular_value), 0.2
        )


class _VariableMassProvider:
    def __init__(self):
        self.base = AnalyticDoubletSOCProviderV220()

    def evaluate_snapshot(self, q):
        base = self.base.evaluate_snapshot(q)
        point = base.point
        mass = np.asarray([[950.0 + 2.0 * float(np.asarray(q)[0])]])
        transformed = ElectronicOperatorPointV21(
            point.q,
            point.H,
            point.dH_dq,
            point.connection_q,
            mass,
            dict(point.metadata),
        ).validate()
        return ElectronicOperatorSnapshotV21(
            transformed, state_vectors=np.eye(4, dtype=complex)
        ).validate()

    def snapshot_overlap(self, left, right):
        return np.eye(4, dtype=complex)


def test_velocity_verlet_rejects_coordinate_dependent_mass():
    with pytest.raises(ValueError, match="requires a constant generalized mass"):
        symmetric_variational_soc_step_v250(_state(), _VariableMassProvider(), 0.2)


def test_static_soc_receipt_without_full_operator_snapshot_is_not_admitted():
    provider = SimpleNamespace(
        evaluate_snapshot=lambda q: SimpleNamespace(matrices=object()),
        snapshot_overlap=lambda left, right: np.eye(4),
    )
    with pytest.raises(TypeError, match="full H, K, D, and mass"):
        symmetric_variational_soc_step_v250(_state(), provider, 0.2)


def test_zero_soc_enabled_and_disabled_trajectories_are_exactly_equivalent():
    enabled = AnalyticSingletTripletSOCProviderV220(
        SingletTripletSOCConfigV220(soc_scale=0.0, soc_enabled=True)
    )
    disabled = AnalyticSingletTripletSOCProviderV220(
        SingletTripletSOCConfigV220(soc_scale=0.0, soc_enabled=False)
    )
    left = run_symmetric_variational_soc_dynamics_v250(
        _state(), enabled, dt_au=0.4, steps=20
    )
    right = run_symmetric_variational_soc_dynamics_v250(
        _state(), disabled, dt_au=0.4, steps=20
    )
    assert np.array_equal(left.final_state.q, right.final_state.q)
    assert np.array_equal(left.final_state.p, right.final_state.p)
    assert np.array_equal(
        left.final_state.electronic_coefficients,
        right.final_state.electronic_coefficients,
    )


def test_step_receipt_rejects_independent_tampering():
    receipt = symmetric_variational_soc_step_v250(
        _state(), AnalyticDoubletSOCProviderV220(), 0.2
    )
    with pytest.raises(ValueError, match="endpoint momentum disagrees"):
        replace(
            receipt,
            end=replace(receipt.end, p=receipt.end.p + 1.0e-3),
        ).validate()
    with pytest.raises(ValueError, match="polar transport disagrees"):
        replace(
            receipt,
            transport_end_to_start=np.diag([1.0, 1.0, 1.0, -1.0]),
        ).validate()
    broken_metrics = dict(receipt.transport_metrics)
    broken_metrics["minimum_singular_value"] = 0.2
    with pytest.raises(ValueError, match="metric minimum_singular_value"):
        replace(receipt, transport_metrics=broken_metrics).validate()
    with pytest.raises(ValueError, match="coordinate-dependent mass"):
        replace(
            receipt,
            mass_matrix_end_au=receipt.mass_matrix_end_au + 1.0e-3,
        ).validate()
