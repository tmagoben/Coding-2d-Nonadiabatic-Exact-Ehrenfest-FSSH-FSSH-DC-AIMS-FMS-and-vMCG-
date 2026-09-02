"""Corrective PySCF static molecular-SOC acceptance campaign for v0.24.1."""

from dataclasses import asdict, dataclass

import numpy as np

from .pyscf_soc_runtime_v241 import run_pyscf_oh_static_soc_evidence_v241
from .pyscf_state_interaction_soc_v241 import (
    BP_SOMF_OPERATOR_FAMILY_V241,
    SpinFreeRootV241,
    assemble_state_interaction_soc_v241,
    clebsch_gordan_twice_v241,
    complete_spin_microstates_v241,
    root_projectors_v241,
    time_reversal_matrix_v241,
)
from .v240_benchmark import run_v0240_release_benchmark


@dataclass(frozen=True)
class V241AcceptanceThresholds:
    expected_inherited_gates: int = 256
    expected_runtime_gates: int = 39
    expected_core_gates: int = 20
    expected_new_gates: int = 59
    expected_total_gates: int = 315


def _raises_v241(callable_object, exceptions, text=None):
    try:
        callable_object()
    except exceptions as exc:
        return text is None or text in str(exc)
    return False


def _core_controls_v241(runtime, inherited):
    doublet_quartet_roots = (
        SpinFreeRootV241("D", -1.0, 1, 1, 0.75),
        SpinFreeRootV241("Q", -0.8, 3, 1, 3.75),
    )
    doublet_quartet_states = complete_spin_microstates_v241(
        doublet_quartet_roots
    )
    J_odd = time_reversal_matrix_v241(doublet_quartet_states)
    odd_identity = np.eye(len(doublet_quartet_states), dtype=complex)

    # A rank-one operator cannot connect S=1/2 and S=5/2.  Populate only that
    # forbidden reduced block and require the assembled direct matrix to remain zero.
    rank_forbidden_roots = (
        SpinFreeRootV241("D", -1.0, 1, 1, 0.75),
        SpinFreeRootV241("X", -0.7, 5, 1, 8.75),
    )
    forbidden_density = np.zeros((2, 2, 2, 2), dtype=complex)
    forbidden_density[0, 1, 0, 1] = 1.0
    forbidden_density[1, 0, 1, 0] = 1.0
    orbital = np.asarray([[0.0, -1j], [1j, 0.0]])
    forbidden_integrals = np.stack((orbital, 0.5 * orbital, -0.25 * orbital))
    forbidden_matrices = assemble_state_interaction_soc_v241(
        rank_forbidden_roots, forbidden_density, forbidden_integrals
    )

    result = runtime.result
    result_payload = result.as_dict(include_large_arrays=False)
    projectors = root_projectors_v241(result.matrices.microstates)
    cross_parity_rejected = _raises_v241(
        lambda: complete_spin_microstates_v241(
            (
                SpinFreeRootV241("S", 0.0, 0, 0),
                SpinFreeRootV241("D", 0.1, 1, 1),
            )
        ),
        (ValueError,),
        "even- and odd-electron",
    )
    wrong_spin_rejected = _raises_v241(
        lambda: SpinFreeRootV241("bad", 0.0, 1, 1, 2.0).validate(),
        (ValueError,),
        "<S^2>",
    )
    noninteger_quantum_number_rejected = _raises_v241(
        lambda: clebsch_gordan_twice_v241(0.5, 0, 2, 0, 2, 0),
        (TypeError,),
        "integer twice",
    )
    return {
        "closed_form_doublet_clebsch_gordan": abs(
            clebsch_gordan_twice_v241(1, 1, 2, 0, 1, 1)
            - 1.0 / np.sqrt(3.0)
        )
        < 1.0e-15,
        "zero_triplet_reference_coefficient_is_exact": (
            clebsch_gordan_twice_v241(2, 0, 2, 0, 2, 0) == 0.0
        ),
        "doublet_quartet_complete_dimension": len(doublet_quartet_states) == 6,
        "doublet_quartet_exact_state_order": tuple(
            state.label for state in doublet_quartet_states
        )
        == (
            "D(M=+1/2)",
            "D(M=-1/2)",
            "Q(M=+3/2)",
            "Q(M=+1/2)",
            "Q(M=-1/2)",
            "Q(M=-3/2)",
        ),
        "odd_time_reversal_is_unitary": np.allclose(
            J_odd.conj().T @ J_odd, odd_identity, atol=1.0e-15
        ),
        "odd_time_reversal_squares_to_minus_identity": np.allclose(
            J_odd @ J_odd.conj(), -odd_identity, atol=1.0e-15
        ),
        "cross_electron_parity_model_space_rejected": cross_parity_rejected,
        "wrong_root_spin_square_rejected": wrong_spin_rejected,
        "noninteger_twice_quantum_number_rejected": (
            noninteger_quantum_number_rejected
        ),
        "rank_one_delta_spin_selection_rule": np.array_equal(
            forbidden_matrices.H_soc,
            np.zeros_like(forbidden_matrices.H_soc),
        ),
        "runtime_evidence_fingerprint_is_sha256": len(runtime.fingerprint()) == 64,
        "provider_result_fingerprint_is_sha256": len(result.fingerprint()) == 64,
        "direct_H_soc_is_serialized_not_reconstructed": (
            "H_soc" in result_payload["matrices"]
            and "soc_eigenvectors" in result_payload["matrices"]
        ),
        "large_ao_tensors_omitted_from_compact_receipt": (
            "effective_ao_cartesian" not in result_payload["integrals"]
        ),
        "root_projectors_cover_direct_matrix": np.array_equal(
            sum(projectors.values(), np.zeros_like(result.matrices.H_soc)),
            np.eye(result.matrices.H_soc.shape[0]),
        ),
        "soc_convention_has_stable_fingerprint": len(
            result.convention.fingerprint()
        )
        == 64,
        "operator_family_is_unambiguously_bp_somf": (
            result.convention.operator_family == BP_SOMF_OPERATOR_FAMILY_V241
        ),
        "static_result_does_not_open_trajectory_admission": (
            not result.molecular_soc_contract.real_backend_admission_ready
        ),
        "prism_is_not_a_runtime_requirement": (
            runtime.claims["Prism_runtime_dependency_required"] is False
        ),
        "inherited_openmolcas_admission_stays_closed": bool(
            inherited["claims"]["external_molecular_SOC_snapshot_admitted"]
            is False
            and inherited["claims"]["native_openmolcas_numeric_crosscheck_implemented"]
            is False
        ),
    }


def run_v0241_release_benchmark(
    thresholds=V241AcceptanceThresholds(),
    *,
    memory_probe_policy="proc_self",
):
    inherited = run_v0240_release_benchmark()
    runtime = run_pyscf_oh_static_soc_evidence_v241(
        memory_probe_policy=memory_probe_policy
    )
    if len(runtime.audit.checks) != thresholds.expected_runtime_gates:
        raise AssertionError("v0.24.1 runtime evidence must define exactly 39 gates.")
    core = _core_controls_v241(runtime, inherited)
    if len(core) != thresholds.expected_core_gates:
        raise AssertionError("v0.24.1 must define exactly 20 core corrective gates.")

    inherited_checks = {
        f"inherited_v0240::{name}": bool(value)
        for name, value in inherited["acceptance"]["checks"].items()
    }
    if len(inherited_checks) != thresholds.expected_inherited_gates:
        raise AssertionError("v0.24.1 must inherit exactly 256 v0.24.0 gates.")
    runtime_checks = {
        f"pyscf_static_soc_runtime::{name}": bool(value)
        for name, value in runtime.audit.checks.items()
    }
    core_checks = {
        f"static_soc_core::{name}": bool(value) for name, value in core.items()
    }
    new_checks = {**runtime_checks, **core_checks}
    if len(new_checks) != thresholds.expected_new_gates:
        raise AssertionError("v0.24.1 must define exactly 59 new gates.")
    checks = {**inherited_checks, **new_checks}
    if len(checks) != thresholds.expected_total_gates:
        raise AssertionError("v0.24.1 must define exactly 315 cumulative gates.")

    return {
        "release": "v0.24.1",
        "theme": (
            "direct static PySCF BP-SOMF molecular SOC with complete spin "
            "microstates and fail-closed trajectory capabilities"
        ),
        "pyscf_static_soc_runtime_evidence": runtime.as_dict(),
        "static_soc_core_controls": core,
        "claims": {
            "real_PySCF_spin_free_runtime_validated": True,
            "real_PySCF_BP_SOMF_execution_validated": True,
            "direct_molecular_SOC_elements_returned": True,
            "doublet_and_Kramers_sector_validated": True,
            "mixed_multiplicity_spin_algebra_validated": True,
            "static_molecular_SOC_tier_validated": True,
            "trajectory_ready_molecular_SOC_validated": False,
            "live_molecular_SOC_backend_admitted": False,
            "physical_SOC_derivatives_validated": False,
            "cross_geometry_SOC_tracking_validated": False,
            "ab_initio_SOC_accuracy_validated": False,
            "external_molecular_SOC_snapshot_admitted": False,
            "native_openmolcas_numeric_crosscheck_implemented": False,
            "Prism_runtime_dependency_required": False,
        },
        "inherited_v0240": inherited,
        "acceptance": {
            "passed": bool(all(checks.values())),
            "checks": checks,
            "inherited_gate_count": len(inherited_checks),
            "runtime_gate_count": len(runtime_checks),
            "core_gate_count": len(core_checks),
            "new_gate_count": len(new_checks),
            "total_gate_count": len(checks),
            "thresholds": asdict(thresholds),
        },
    }
