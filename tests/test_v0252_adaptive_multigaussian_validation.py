import pytest

from gaussian_dynamics.adaptive_multigaussian_tdvp_validation_v252 import (
    V252_CONVERGENCE_DT_AU,
    run_adaptive_multigaussian_validation_evidence_v252,
)


@pytest.fixture(scope="module")
def evidence():
    return run_adaptive_multigaussian_validation_evidence_v252()


def test_validation_evidence_passes_70_native_boolean_gates(evidence):
    assert evidence.audit.passed
    assert len(evidence.audit.checks) == 70
    assert all(type(value) is bool for value in evidence.audit.checks.values())
    assert all(evidence.audit.checks.values())


def test_width_chirp_evolution_reversal_and_covariance_are_resolved(evidence):
    metrics = evidence.audit.metrics
    assert metrics["even_width_change"] > 1.0e-5
    assert metrics["odd_chirp_change"] > 1.0e-4
    assert metrics["even_reversibility"]["widths"] < 2.0e-13
    assert metrics["odd_reversibility"]["chirps"] < 2.0e-13
    assert metrics["permutation_errors"]["widths"] < 1.0e-14
    assert metrics["constant_gauge_errors"]["coefficients"] < 1.0e-13


def test_duplicate_null_space_and_harmonic_oracles_are_explicit(evidence):
    null_receipt = evidence.audit.metrics["compatible_null_space"]
    assert null_receipt["rank"] == 6
    assert null_receipt["nullity"] == 6
    assert null_receipt["null_rhs_relative"] < 1.0e-14
    assert null_receipt["linear_residual_relative"] < 2.0e-14
    harmonic = evidence.audit.metrics["harmonic_errors"]
    assert harmonic["q"] < 2.0e-8
    assert harmonic["p"] < 2.0e-8
    assert harmonic["width"] < 2.0e-8
    assert harmonic["chirp"] < 6.0e-8


def test_adaptive_implicit_midpoint_refinement_is_second_order(evidence):
    assert V252_CONVERGENCE_DT_AU == (0.1, 0.05, 0.025, 0.0125)
    ratios = evidence.audit.metrics["convergence_state_change_ratios"]
    assert len(ratios) == 2
    assert all(0.245 < value < 0.255 for value in ratios)


def test_evidence_fingerprint_is_deterministic(evidence):
    assert len(evidence.fingerprint()) == 64
    assert evidence.fingerprint() == evidence.fingerprint()

