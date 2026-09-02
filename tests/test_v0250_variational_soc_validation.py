import pytest

from gaussian_dynamics.variational_soc_validation_v250 import (
    V250_CONVERGENCE_DT_AU,
    run_variational_soc_validation_evidence_v250,
)


@pytest.fixture(scope="module")
def evidence():
    return run_variational_soc_validation_evidence_v250()


def test_validation_evidence_passes_45_native_boolean_gates(evidence):
    assert evidence.audit.passed
    assert len(evidence.audit.checks) == 45
    assert all(type(value) is bool for value in evidence.audit.checks.values())
    assert all(evidence.audit.checks.values())


def test_even_odd_and_complex_gauge_metrics_are_machine_scale(evidence):
    metrics = evidence.audit.metrics
    assert metrics["even_maximum_norm_drift"] < 1.0e-13
    assert metrics["odd_maximum_norm_drift"] < 1.0e-13
    assert metrics["even_reversibility"]["spinor"] < 1.0e-13
    assert metrics["odd_reversibility"]["spinor"] < 1.0e-13
    assert metrics["gauge_errors"]["spinor"] < 1.0e-13


def test_timestep_and_energy_errors_share_second_order_plateau(evidence):
    assert V250_CONVERGENCE_DT_AU == (0.8, 0.4, 0.2, 0.1)
    state_ratios = evidence.audit.metrics["convergence_state_change_ratios"]
    energy_ratios = evidence.audit.metrics["convergence_energy_drift_ratios"]
    assert all(0.249 < value < 0.251 for value in state_ratios)
    assert all(0.249 < value < 0.251 for value in energy_ratios)


def test_raw_contraction_and_unitary_polar_transport_remain_distinct(evidence):
    metrics = evidence.audit.metrics
    assert metrics["minimum_retained_singular_value"] == pytest.approx(0.97)
    assert metrics["maximum_transport_unitarity_residual"] < 2.0e-14


def test_full_tdvp_and_real_pyscf_trajectory_claims_remain_closed(evidence):
    claims = evidence.claims
    assert claims["restricted_single_packet_tdvp_validated"] is True
    assert claims["symmetric_strang_verlet_coupling_validated"] is True
    assert claims["svd_computed_polar_transport_validated"] is True
    assert claims["full_multi_gaussian_tdvp_validated"] is False
    assert claims["adaptive_gaussian_width_tdvp_validated"] is False
    assert claims["plain_verlet_for_general_tdvp_validated"] is False
    assert claims["coordinate_dependent_mass_verlet_validated"] is False
    assert claims["real_pyscf_soc_trajectory_admitted"] is False
    assert claims["general_ab_initio_soc_dynamics_accuracy_validated"] is False
