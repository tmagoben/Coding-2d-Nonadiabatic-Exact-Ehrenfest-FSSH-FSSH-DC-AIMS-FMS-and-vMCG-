import os

import pytest

from gaussian_dynamics.pyscf_runtime_v232 import (
    probe_pyscf_runtime_v232,
    run_pyscf_runtime_evidence_v232,
)


_PROBE = probe_pyscf_runtime_v232()
pytestmark = pytest.mark.skipif(
    not _PROBE.usable,
    reason="real PySCF 2.13.1 runtime integration test",
)


@pytest.fixture(scope="module")
def evidence():
    import pyscf

    original_current_memory = pyscf.lib.current_memory
    result = run_pyscf_runtime_evidence_v232()
    assert pyscf.lib.current_memory is original_current_memory
    return result


def test_real_runtime_provenance_and_memory_workaround(evidence):
    assert evidence.passed
    assert evidence.runtime.pyscf_distribution_version == "2.13.1"
    assert evidence.runtime.pyscf_module_version == "2.13.1"
    assert evidence.runtime.pyscf_distribution_file_count > 1000
    assert evidence.runtime.pyscf_verified_file_count > 1000
    assert evidence.runtime.pyscf_verified_size_bytes > 100_000_000
    assert len(evidence.runtime.pyscf_verified_content_sha256) == 64
    assert len(evidence.runtime.environment_sha256) == 64

    pid_statm_exists = os.path.exists(f"/proc/{os.getpid()}/statm")
    if not pid_statm_exists:
        assert evidence.runtime.memory_probe_mode == (
            "proc_self_statm_pid_namespace_fallback"
        )


def test_real_two_geometry_energy_gradient_nac_and_overlap_evidence(evidence):
    metrics = evidence.smoke["metrics"]

    assert evidence.checks["two_geometry_values_are_finite"]
    assert evidence.checks["analytic_gradients_are_nontrivial"]
    assert evidence.checks["nac_is_nontrivial"]
    assert evidence.checks["nac_is_antisymmetric"]
    assert evidence.checks["self_overlap_is_identity"]
    assert evidence.checks["cross_overlap_is_reciprocal"]
    assert evidence.checks["cross_overlap_is_a_physical_contraction"]
    assert evidence.checks["cross_overlap_is_not_forced_to_exact_isometry"]
    assert metrics["minimum_singular_value"] > 0.99
    assert metrics["cross_overlap_isometry_defect"] > 1e-8


def test_real_nac_mapping_matches_phase_aligned_central_differences(evidence):
    metrics = evidence.nac_mapping["metrics"]
    production_errors = metrics["production_central_difference_max_errors"]
    direct_errors = metrics["direct_state_i_j_central_difference_errors"]

    assert evidence.checks["production_state_tuple_mapping_has_correct_sign"]
    assert evidence.checks["production_mapping_shows_second_order_convergence"]
    assert evidence.checks["raw_state_i_j_matches_overlap_derivative"]
    assert evidence.checks[
        "raw_state_i_j_mapping_shows_second_order_convergence"
    ]
    assert production_errors[2] < production_errors[1] < production_errors[0]
    assert direct_errors[2] < direct_errors[1] < direct_errors[0]
    # The absolute root phase is arbitrary and can flip with BLAS threading.
    # The certified statement is the sign *relative to the phase-aligned
    # overlap derivative*, not a hard-coded positive/negative number.
    assert (
        metrics["production_dij_selected"]
        * metrics["finest_production_fd_selected_real"]
        > 0.0
    )
    assert (
        metrics["raw_state_i_j_selected"]
        * metrics["raw_state_j_i_selected"]
        < 0.0
    )


def test_real_etf_and_scaled_nac_semantics_are_separate(evidence):
    assert evidence.checks["opposite_state_tuples_are_antisymmetric_no_etf"]
    assert evidence.checks["opposite_state_tuples_are_antisymmetric_etf"]
    assert evidence.checks["scaled_nac_is_symmetric_under_tuple_swap"]
    assert evidence.checks["scaled_nac_obeys_energy_difference_relation"]
    assert evidence.checks["production_scaled_nac_is_symmetric"]
    assert evidence.checks[
        "production_scaled_nac_obeys_energy_difference_relation"
    ]
    assert evidence.checks["full_overlap_derivative_uses_no_etf"]
    assert evidence.checks["etf_removes_translation_component"]
    assert evidence.checks["etf_and_full_overlap_derivatives_are_distinct"]
