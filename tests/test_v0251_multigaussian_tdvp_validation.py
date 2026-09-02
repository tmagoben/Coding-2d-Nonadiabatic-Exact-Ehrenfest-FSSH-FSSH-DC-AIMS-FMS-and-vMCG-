import pytest

from gaussian_dynamics.multigaussian_tdvp_validation_v251 import (
    V251_CONVERGENCE_DT_AU,
    run_multigaussian_tdvp_validation_evidence_v251,
)


@pytest.fixture(scope="module")
def evidence():
    return run_multigaussian_tdvp_validation_evidence_v251()


def test_validation_evidence_passes_55_native_boolean_gates(evidence):
    assert evidence.audit.passed
    assert len(evidence.audit.checks) == 55
    assert all(type(value) is bool for value in evidence.audit.checks.values())
    assert all(evidence.audit.checks.values())


def test_even_odd_reversal_and_covariance_errors_are_machine_scale(evidence):
    metrics = evidence.audit.metrics
    assert metrics["even_reversibility"]["coefficients"] < 1.0e-13
    assert metrics["odd_reversibility"]["coefficients"] < 1.0e-13
    assert metrics["permutation_errors"]["coefficients"] < 1.0e-13
    assert metrics["constant_gauge_errors"]["coefficients"] < 1.0e-13
    assert metrics["even_maximum_norm_drift"] < 1.0e-12
    assert metrics["odd_maximum_norm_drift"] < 1.0e-12


def test_compatible_null_space_and_harmonic_reduction_are_explicit(evidence):
    null_receipt = evidence.audit.metrics["compatible_null_space"]
    assert null_receipt["rank"] == 4
    assert null_receipt["nullity"] == 4
    assert null_receipt["null_rhs_relative"] < 1.0e-14
    assert null_receipt["linear_residual_relative"] < 2.0e-14
    assert evidence.audit.metrics["harmonic_reduction_errors"]["qdot"] < 1.0e-14
    assert evidence.audit.metrics["harmonic_reduction_errors"]["pdot"] < 1.0e-14


def test_implicit_midpoint_refinement_is_second_order(evidence):
    assert V251_CONVERGENCE_DT_AU == (0.1, 0.05, 0.025, 0.0125)
    ratios = evidence.audit.metrics["convergence_state_change_ratios"]
    assert len(ratios) == 2
    assert all(0.245 < value < 0.255 for value in ratios)


def test_evidence_fingerprint_is_deterministic(evidence):
    assert len(evidence.fingerprint()) == 64
    assert evidence.fingerprint() == evidence.fingerprint()

