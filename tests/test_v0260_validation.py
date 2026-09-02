import pytest

from gaussian_dynamics.multidimensional_validation_v260 import (
    MULTIDIMENSIONAL_VALIDATION_SCHEMA_V260,
    run_multidimensional_validation_evidence_v260,
)


@pytest.fixture(scope="module")
def evidence():
    return run_multidimensional_validation_evidence_v260()


def test_v0260_validation_passes_all_80_gates(evidence):
    assert evidence.passed is True
    assert evidence.check_count == 80
    assert len(evidence.checks) == 80
    assert all(evidence.checks.values())


def test_v0260_validation_schema_and_fingerprints_are_stable(evidence):
    assert MULTIDIMENSIONAL_VALIDATION_SCHEMA_V260.endswith("v0.26.0")
    assert len(evidence.exact_grid_fingerprint) == 64
    assert len(evidence.controlled_trajectory_fingerprint) == 64
    assert len(evidence.fingerprint()) == 64


def test_v0260_adaptation_improves_exact_reference_error(evidence):
    metrics = evidence.metrics
    assert metrics["controlled_exact_wavefunction_error"] < metrics["one_packet_exact_wavefunction_error"]
    assert metrics["controlled_density_error"] < metrics["one_packet_density_error"]


def test_v0260_claim_boundaries_remain_honest(evidence):
    claims = evidence.as_dict()["claims"]
    assert claims["tdvp"]["full_correlated_width_matrices_validated"] is False
    assert claims["tdvp"]["coordinate_dependent_electronic_gauge_covariance_validated"] is False
    assert claims["basis"]["full_aims_branching_validated"] is False
    assert claims["basis"]["real_pyscf_soc_trajectory_admitted"] is False
    assert claims["exact_grid"]["absorbing_boundary_conditions_validated"] is False
