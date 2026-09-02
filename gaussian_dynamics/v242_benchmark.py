"""Connected-geometry PySCF SOC differential acceptance campaign for v0.24.2."""

from dataclasses import asdict, dataclass, replace
import hashlib
import inspect
import json

import numpy as np

from .finite_manifold_transport_v233 import certified_transport_from_overlap_v233
from .pyscf_differential_soc_v242 import (
    OH_BOND_STEPS_BOHR_V242,
    PYSCF_DIFFERENTIAL_SOC_CAPABILITY_V242,
    audit_pyscf_oh_bond_differential_soc_v242,
    build_pyscf_bp_somf_integrals_direct_jk_v242,
    phase_align_complete_multiplet_overlap_v242,
    run_pyscf_oh_bond_differential_evidence_v242,
)
from .v241_benchmark import run_v0241_release_benchmark


@dataclass(frozen=True)
class V242AcceptanceThresholds:
    expected_inherited_gates: int = 315
    expected_runtime_gates: int = 60
    expected_core_gates: int = 25
    expected_new_gates: int = 85
    expected_total_gates: int = 400


def _raises_v242(callable_object, exceptions, text=None):
    try:
        callable_object()
    except exceptions as exc:
        return text is None or text in str(exc)
    return False


def _core_controls_v242(evidence):
    scan = evidence.scan
    center = scan.center
    records = scan.derivative_records
    claims = evidence.claims
    microstates = center.matrices.microstates
    source = inspect.getsource(build_pyscf_bp_somf_integrals_direct_jk_v242)

    phases = np.asarray(
        [1j, 1j, -1.0, -1.0, np.exp(0.37j), np.exp(0.37j)],
        dtype=complex,
    )
    aligned, corrections = phase_align_complete_multiplet_overlap_v242(
        np.diag(phases), microstates
    )
    degenerate_swap = np.kron(
        np.asarray(
            [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
            dtype=complex,
        ),
        np.eye(2, dtype=complex),
    )

    tampered_direct_jk = replace(
        scan,
        direct_jk_explicit_max_abs_error=1.0e-4,
    )
    tampered_direct_jk_audit = audit_pyscf_oh_bond_differential_soc_v242(
        tampered_direct_jk
    )
    tampered_claims = dict(scan.claims)
    tampered_claims["full_cartesian_derivative_tensor_validated"] = True
    tampered_claim_audit = audit_pyscf_oh_bond_differential_soc_v242(
        replace(scan, claims=tampered_claims)
    )

    malformed_soc = records[0].K_soc.copy()
    malformed_soc[0, 1] += 1.0e-3
    compensated_shift = np.zeros_like(records[0].K_total)
    compensated_shift[0, 0] = 1.0e-3
    nonhermitian_rejected = _raises_v242(
        lambda: replace(records[0], K_soc=malformed_soc).validate(),
        (ValueError,),
        "K_soc is not Hermitian",
    )
    compensated_tamper_rejected = _raises_v242(
        lambda: replace(
            records[0],
            K_spin_free=records[0].K_spin_free + compensated_shift,
            K_soc=records[0].K_soc - compensated_shift,
        ).validate(),
        (ValueError,),
        "spin-free derivative disagrees",
    )

    canonical_payload = json.dumps(
        evidence.as_dict(),
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
    ).encode("utf-8")
    canonical_digest = hashlib.sha256(canonical_payload).hexdigest()
    J = np.asarray(center.matrices.time_reversal_matrix, dtype=complex)
    identity = np.eye(J.shape[0], dtype=complex)

    return {
        "production_source_calls_direct_jk": "jk.get_jk" in source,
        "production_source_avoids_rank_five_mol_intor": (
            "mol.intor(BP_SOMF_TWO_ELECTRON_INTEGRAL_V241" not in source
        ),
        "canonical_second_order_step_ladder_frozen": (
            tuple(record.displacement_bohr for record in records)
            == OH_BOND_STEPS_BOHR_V242
        ),
        "root_phase_is_shared_across_each_doublet": bool(
            np.allclose(aligned, identity)
            and np.allclose(corrections[0], corrections[1])
            and np.allclose(corrections[2], corrections[3])
            and np.allclose(corrections[4], corrections[5])
        ),
        "degenerate_rotation_is_rejected_by_phase_only_rule": _raises_v242(
            lambda: phase_align_complete_multiplet_overlap_v242(
                degenerate_swap, microstates
            ),
            (ValueError,),
            "phase is ambiguous",
        ),
        "three_doublets_lift_to_six_microstates": len(microstates) == 6,
        "complete_doublet_time_reversal_squares_to_minus_identity": bool(
            np.allclose(J @ J.conj(), -identity, atol=1.0e-14)
        ),
        "spectral_expansion_overlap_is_rejected": _raises_v242(
            lambda: certified_transport_from_overlap_v233(
                np.diag([1.01, 0.9]).astype(complex)
            ),
            (ValueError,),
            "physically inconsistent",
        ),
        "rank_lost_overlap_is_rejected": _raises_v242(
            lambda: certified_transport_from_overlap_v233(
                np.diag([1.0, 0.0]).astype(complex)
            ),
            (ValueError,),
            "not trajectory ready",
        ),
        "tampered_direct_jk_oracle_error_fails_audit": bool(
            not tampered_direct_jk_audit.passed
            and not tampered_direct_jk_audit.checks[
                "direct_jk_matches_explicit_tensor_oracle"
            ]
        ),
        "inflated_full_cartesian_claim_fails_audit": bool(
            not tampered_claim_audit.passed
            and not tampered_claim_audit.checks[
                "full_cartesian_derivative_claim_remains_false"
            ]
        ),
        "nonhermitian_derivative_record_is_rejected": nonhermitian_rejected,
        "compensated_component_tamper_is_endpoint_rejected": (
            compensated_tamper_rejected
        ),
        "nondecreasing_step_sequence_is_rejected": _raises_v242(
            lambda: replace(
                scan,
                derivative_records=(records[0], records[1], records[1]),
            ).validate(),
            (ValueError,),
            "strictly decreasing",
        ),
        "evidence_fingerprint_is_sha256": bool(
            len(evidence.fingerprint()) == 64
            and all(character in "0123456789abcdef" for character in evidence.fingerprint())
        ),
        "canonical_evidence_serialization_is_stable": bool(
            evidence.fingerprint() == canonical_digest
            and canonical_digest
            == hashlib.sha256(canonical_payload).hexdigest()
        ),
        "capability_is_connected_geometry_preview": bool(
            PYSCF_DIFFERENTIAL_SOC_CAPABILITY_V242
            == "connected_geometry_differential_preview"
        ),
        "direct_jk_execution_claim_is_true": bool(
            claims["direct_jk_somf_execution_validated"] is True
        ),
        "continuous_physical_connection_claim_is_false": bool(
            claims["continuous_physical_derivative_connection_validated"] is False
        ),
        "full_cartesian_derivative_claim_is_false": bool(
            claims["full_cartesian_derivative_tensor_validated"] is False
        ),
        "analytic_soc_derivative_claim_is_false": bool(
            claims["analytic_soc_derivatives_validated"] is False
        ),
        "real_mixed_multiplicity_runtime_claim_is_false": bool(
            claims["real_mixed_multiplicity_runtime_validated"] is False
        ),
        "trajectory_ready_claim_is_false": bool(
            claims["trajectory_ready_molecular_soc_validated"] is False
        ),
        "live_backend_admission_claim_is_false": bool(
            claims["live_molecular_soc_backend_admitted"] is False
        ),
        "ab_initio_accuracy_claim_is_false": bool(
            claims["ab_initio_soc_accuracy_validated"] is False
        ),
    }


def run_v0242_release_benchmark(
    thresholds=V242AcceptanceThresholds(),
    *,
    memory_probe_policy="proc_self",
):
    inherited = run_v0241_release_benchmark(
        memory_probe_policy=memory_probe_policy
    )
    evidence = run_pyscf_oh_bond_differential_evidence_v242(
        memory_probe_policy=memory_probe_policy
    )
    if len(evidence.audit.checks) != thresholds.expected_runtime_gates:
        raise AssertionError("v0.24.2 runtime evidence must define exactly 60 gates.")
    core = _core_controls_v242(evidence)
    if len(core) != thresholds.expected_core_gates:
        raise AssertionError("v0.24.2 must define exactly 25 core gates.")

    inherited_checks = {
        f"inherited_v0241::{name}": bool(value)
        for name, value in inherited["acceptance"]["checks"].items()
    }
    if len(inherited_checks) != thresholds.expected_inherited_gates:
        raise AssertionError("v0.24.2 must inherit exactly 315 v0.24.1 gates.")
    runtime_checks = {
        f"pyscf_differential_soc_runtime::{name}": bool(value)
        for name, value in evidence.audit.checks.items()
    }
    core_checks = {
        f"differential_soc_core::{name}": bool(value)
        for name, value in core.items()
    }
    new_checks = {**runtime_checks, **core_checks}
    if len(new_checks) != thresholds.expected_new_gates:
        raise AssertionError("v0.24.2 must define exactly 85 new gates.")
    checks = {**inherited_checks, **new_checks}
    if len(checks) != thresholds.expected_total_gates:
        raise AssertionError("v0.24.2 must define exactly 400 cumulative gates.")

    return {
        "release": "v0.24.2",
        "theme": (
            "direct-JK BP-SOMF, connected PySCF geometry snapshots, complete-"
            "multiplet polar transport, and finite-difference SOC components"
        ),
        "pyscf_differential_soc_evidence": evidence.as_dict(),
        "differential_soc_core_controls": core,
        "claims": dict(evidence.claims),
        "inherited_v0241": inherited,
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
