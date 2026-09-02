import pytest

from gaussian_dynamics.correlated_validation_v270 import (
    CORRELATED_VALIDATION_SCHEMA_V270,
    run_correlated_validation_evidence_v270,
)


@pytest.fixture(scope="module")
def evidence():
    return run_correlated_validation_evidence_v270()


def test_v0270_validation_passes_all_100_gates(evidence):
    assert evidence.passed is True
    assert evidence.check_count == 100
    assert len(evidence.checks) == 100
    assert all(evidence.checks.values())


def test_v0270_validation_schema_and_fingerprints_are_stable(evidence):
    assert CORRELATED_VALIDATION_SCHEMA_V270.endswith("v0.27.0")
    assert len(evidence.trajectory_fingerprint) == 64
    assert len(evidence.lifecycle_fingerprint) == 64
    assert len(evidence.fingerprint()) == 64


def test_v0270_exact_riccati_oracle_shows_second_order_convergence(evidence):
    metrics = evidence.metrics
    assert metrics["riccati_error_dt_002"] < metrics["riccati_error_dt_004"]
    assert metrics["riccati_error_dt_001"] < metrics["riccati_error_dt_002"]
    assert metrics["riccati_order_coarse"] > 1.95
    assert metrics["riccati_order_fine"] > 1.95


def test_v0270_claim_boundaries_remain_honest(evidence):
    claims = evidence.as_dict()["claims"]
    assert claims["tdvp"]["full_correlated_width_matrices_validated"] is True
    assert claims["tdvp"]["arbitrary_orthogonal_coordinate_covariance_validated"] is True
    assert claims["tdvp"]["coordinate_dependent_electronic_gauge_covariance_validated"] is False
    assert claims["tdvp"]["real_pyscf_soc_trajectory_admitted"] is False
    assert claims["basis"]["full_aims_branching_validated"] is False
