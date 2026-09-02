import inspect
from dataclasses import replace

import numpy as np
import pytest

from gaussian_dynamics.pyscf_differential_soc_v242 import (
    OH_BOND_STEPS_BOHR_V242,
    PYSCF_DIRECT_JK_SOMF_STRATEGY_V242,
    audit_pyscf_oh_bond_differential_soc_v242,
    build_pyscf_bp_somf_integrals_direct_jk_v242,
    run_pyscf_oh_bond_differential_evidence_v242,
)
from gaussian_dynamics.pyscf_state_interaction_soc_v241 import (
    probe_pyscf_static_soc_runtime_v241,
)


@pytest.fixture(scope="module")
def evidence():
    probe = probe_pyscf_static_soc_runtime_v241()
    if not probe.usable:
        pytest.skip(probe.failure_reason)
    return run_pyscf_oh_bond_differential_evidence_v242()


def test_direct_jk_runtime_path_matches_explicit_tensor_oracle(evidence):
    assert evidence.audit.passed
    assert len(evidence.audit.checks) == 60
    assert all(evidence.audit.checks.values())
    assert evidence.scan.direct_jk_explicit_max_abs_error < 2.0e-14
    assert (
        evidence.scan.center.calculation_input["somf_contraction"]
        == PYSCF_DIRECT_JK_SOMF_STRATEGY_V242
    )
    source = inspect.getsource(build_pyscf_bp_somf_integrals_direct_jk_v242)
    assert "mol.intor(BP_SOMF_TWO_ELECTRON_INTEGRAL_V241" not in source
    assert "jk.get_jk" in source


def test_oh_connected_geometry_overlaps_are_contractive_and_well_retained(evidence):
    records = evidence.scan.derivative_records
    assert len(evidence.scan.endpoint_snapshots) == 3
    assert sum(len(pair) for pair in evidence.scan.endpoint_snapshots) == 6
    assert tuple(record.displacement_bohr for record in records) == OH_BOND_STEPS_BOHR_V242
    blocks = [
        record.overlap_metrics[side]
        for record in records
        for side in ("minus", "plus")
    ]
    assert min(block["minimum_singular_value"] for block in blocks) > 0.996
    assert max(block["maximum_singular_value"] for block in blocks) <= 1.0 + 1.0e-10
    assert all(block["physically_consistent"] for block in blocks)
    assert all(block["trajectory_ready"] for block in blocks)


def test_endpoint_receipts_are_serialized_and_fingerprint_bound(evidence):
    payload = evidence.as_dict()["scan"]
    assert len(payload["endpoint_snapshots"]) == 3
    assert set(payload["endpoint_snapshots"][0]) == {"minus", "plus"}

    first = evidence.scan.derivative_records[0]
    tampered = replace(first, minus_fingerprint="0" * 64)
    broken_scan = replace(
        evidence.scan,
        derivative_records=(tampered,) + evidence.scan.derivative_records[1:],
    )
    with pytest.raises(ValueError, match="minus endpoint fingerprint is not bound"):
        broken_scan.validate()


def test_component_derivatives_are_separate_hermitian_and_time_reversal_safe(evidence):
    for record in evidence.scan.derivative_records:
        assert np.allclose(
            record.K_total, record.K_spin_free + record.K_soc, atol=1.0e-12
        )
        assert np.allclose(record.K_spin_free, record.K_spin_free.conj().T)
        assert np.allclose(record.K_soc, record.K_soc.conj().T)
        assert record.residuals["K_soc_time_reversal"] < 1.0e-12
        assert record.residuals["D_antihermiticity"] == 0.0


def test_spin_free_and_soc_derivatives_share_second_order_plateau(evidence):
    metrics = evidence.scan.convergence_metrics
    assert 0.24 < metrics["K_spin_free"]["fine_to_coarse_change_ratio"] < 0.27
    assert 0.24 < metrics["K_soc"]["fine_to_coarse_change_ratio"] < 0.28
    assert metrics["K_soc"]["finest_norm_frobenius"] > 1.0e-4
    slopes = [
        record.residuals["parallel_transport_hermitian_slope_frobenius"]
        for record in evidence.scan.derivative_records
    ]
    assert slopes[0] > slopes[1] > slopes[2]


def test_differential_preview_keeps_unvalidated_claims_closed(evidence):
    claims = evidence.claims
    assert claims["direct_jk_somf_execution_validated"] is True
    assert claims["connected_geometry_soc_snapshots_validated"] is True
    assert claims["transported_soc_derivative_preview_validated"] is True
    assert claims["continuous_physical_derivative_connection_validated"] is False
    assert claims["full_cartesian_derivative_tensor_validated"] is False
    assert claims["analytic_soc_derivatives_validated"] is False
    assert claims["real_mixed_multiplicity_runtime_validated"] is False
    assert claims["trajectory_ready_molecular_soc_validated"] is False
    assert claims["live_molecular_soc_backend_admitted"] is False
    assert claims["ab_initio_soc_accuracy_validated"] is False


def test_audit_recomputes_same_60_native_boolean_gates(evidence):
    repeated = audit_pyscf_oh_bond_differential_soc_v242(evidence.scan)
    assert repeated.checks == evidence.audit.checks
    assert all(type(value) is bool for value in repeated.checks.values())
