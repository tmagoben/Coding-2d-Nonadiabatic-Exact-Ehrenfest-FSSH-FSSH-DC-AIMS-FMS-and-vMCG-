"""Real PySCF runtime and admission-correction campaign for v0.23.2."""

from dataclasses import asdict, dataclass, replace
from pathlib import Path
from types import SimpleNamespace
import tempfile

import numpy as np

from .molecular_soc_admission_v232 import audit_molecular_soc_provider_v232
from .molecular_soc_runtime_v232 import (
    BackendAdmissionPolicyV232,
    BackendMethodIdentityV232,
    CONVERGENCE_METADATA_KEY_V232,
    ConvergenceMetadataV232,
    convergence_from_snapshot_v232,
)
from .pyscf_nac_convention_v232 import (
    PYSCF_NAC_EMPIRICAL_MAPPING_V232,
    PYSCF_NAC_UPSTREAM_DOCUMENTATION_V232,
    pyscf_state_tuple_for_internal_dij_v232,
)
from .pyscf_runtime_v232 import run_pyscf_runtime_evidence_v232
from .pyscf_soc_adapter_v231 import require_pyscf_soc_adapter_v231
from .pyscf_soc_adapter_v232 import (
    PYSCF_NAC_CONVENTION_V232,
    PySCFMethodSpecificCapabilitiesV232,
    validate_pyscf_engine_contract_v232,
)
from .v230_benchmark import build_v230_reference_replay
from .v231_benchmark import build_v231_admission_bundle, run_v0231_release_benchmark


@dataclass(frozen=True)
class V232AcceptanceThresholds:
    expected_inherited_gates: int = 123
    expected_runtime_gates: int = 28
    expected_new_gates: int = 45
    expected_total_gates: int = 168


def _raises_v232(callable_object, exceptions, text=None):
    try:
        callable_object()
    except exceptions as exc:
        return text is None or text in str(exc)
    return False


def _overlap_controls_v232(directory):
    dataset = build_v230_reference_replay(Path(directory) / "reference")
    overlaps = np.asarray(dataset.overlaps, dtype=complex).copy()
    nrecord, _, nstate, _ = overlaps.shape
    contraction = 0.9 * np.eye(nstate, dtype=complex)
    for left in range(nrecord):
        for right in range(left + 1, nrecord):
            overlaps[left, right] = contraction
            overlaps[right, left] = contraction.conj().T
    contracted = replace(dataset, overlaps=overlaps)
    contracted_accepted = contracted.validate() is contracted
    diagnostics = contracted.overlap_diagnostics()
    diagnostics_physical = bool(
        diagnostics.maximum_self_identity_residual == 0.0
        and diagnostics.maximum_reciprocity_residual == 0.0
        and abs(diagnostics.minimum_cross_geometry_singular_value - 0.9) < 1e-12
        and abs(diagnostics.maximum_cross_geometry_singular_value - 0.9) < 1e-12
        and diagnostics.maximum_contraction_excess == 0.0
    )

    expansive_table = overlaps.copy()
    expansive = np.zeros((nstate, nstate), dtype=complex)
    expansive[0, 0] = 0.8
    expansive[0, 1] = 0.8
    expansive_table[0, 1] = expansive
    expansive_table[1, 0] = expansive.conj().T
    expansion_rejected = _raises_v232(
        lambda: replace(dataset, overlaps=expansive_table).validate(),
        (ValueError,),
        "singular values",
    )

    nonreciprocal = overlaps.copy()
    nonreciprocal[0, 1] = 0.9 * np.eye(nstate)
    nonreciprocal[1, 0] = 0.8 * np.eye(nstate)
    reciprocity_rejected = _raises_v232(
        lambda: replace(dataset, overlaps=nonreciprocal).validate(),
        (ValueError,),
        "reciprocity",
    )

    bad_self = overlaps.copy()
    bad_self[0, 0] *= 0.999
    self_identity_rejected = _raises_v232(
        lambda: replace(dataset, overlaps=bad_self).validate(),
        (ValueError,),
        "identity",
    )
    return {
        "physical_contraction_accepted": contracted_accepted,
        "physical_contraction_diagnostics": diagnostics_physical,
        "spectral_expansion_rejected": expansion_rejected,
        "reciprocity_corruption_rejected": reciprocity_rejected,
        "self_identity_corruption_rejected": self_identity_rejected,
        "diagnostics": diagnostics.as_dict(),
    }


def _trusted_pyscf_identity_v232(**changes):
    fields = {
        "backend_name": "PySCF",
        "backend_version": "2.13.1",
        "source_kind": "live_ab_initio",
        "adapter_name": "v232-negative-control-adapter",
        "adapter_version": "1",
        "electronic_method": "SA-CASSCF/SI-SOC negative-control identity",
        "basis": "test basis",
        "active_space": "CAS(2,2)",
        "soc_operator": "test SOC operator",
        "scalar_relativistic_method": "none",
        "derivative_method": "test derivatives",
        "nac_convention": PYSCF_NAC_CONVENTION_V232,
    }
    fields.update(changes)
    return BackendMethodIdentityV232(**fields).validate()


class _CompleteEngineV232:
    def __init__(self, identity, *, complete=True, noncallable=False):
        self.method_identity = identity
        self.capabilities = PySCFMethodSpecificCapabilitiesV232(
            state_interaction_soc=bool(complete),
            physical_soc_derivatives=True,
            analytic_spin_free_gradients=True,
            derivative_connections=True,
            many_electron_overlaps=True,
            raw_artifact_parser=True,
            fresh_execution=True,
        )
        if noncallable:
            self.components = None

    def components(self, q):
        raise AssertionError("negative-control engine must not execute")

    def evaluate_snapshot(self, q):
        raise AssertionError("negative-control engine must not execute")

    def snapshot_overlap(self, left, right):
        raise AssertionError("negative-control engine must not execute")

    def write_raw_artifacts(self, *args, **kwargs):
        raise AssertionError("negative-control engine must not execute")

    def validate_raw_artifacts_v232(self, *args, **kwargs):
        raise AssertionError("negative-control engine must not execute")


def _admission_controls_v232(directory):
    identity = _trusted_pyscf_identity_v232()

    class DeclarationTrap:
        method_identity = identity

        @property
        def capabilities(self):
            raise AssertionError("capabilities consulted before structure")

    structure_before_flags = _raises_v232(
        lambda: validate_pyscf_engine_contract_v232(
            DeclarationTrap(),
            runtime_version="2.13.1",
            expected_identity=identity,
        ),
        (TypeError,),
        "incomplete",
    )
    noncallable_rejected = _raises_v232(
        lambda: validate_pyscf_engine_contract_v232(
            _CompleteEngineV232(identity, noncallable=True),
            runtime_version="2.13.1",
            expected_identity=identity,
        ),
        (TypeError,),
        "noncallable",
    )
    identity_mismatch_rejected = _raises_v232(
        lambda: validate_pyscf_engine_contract_v232(
            _CompleteEngineV232(identity),
            runtime_version="2.13.1",
            expected_identity=_trusted_pyscf_identity_v232(basis="other basis"),
        ),
        (ValueError,),
        "trusted identity",
    )
    incomplete_capabilities_rejected = _raises_v232(
        lambda: validate_pyscf_engine_contract_v232(
            _CompleteEngineV232(identity, complete=False),
            runtime_version="2.13.1",
            expected_identity=identity,
        ),
        (RuntimeError,),
        "incomplete",
    )

    legacy_snapshot = SimpleNamespace(
        point=SimpleNamespace(metadata={"scf_converged": True}),
        metadata={"derivatives_converged": True},
    )
    legacy_convergence_rejected = _raises_v232(
        lambda: convergence_from_snapshot_v232(legacy_snapshot),
        (ValueError,),
        "lacks canonical",
    )
    convergence = ConvergenceMetadataV232(
        scf=True,
        correlated_wavefunction=True,
        state_interaction_soc=True,
        spin_free_gradients=True,
        soc_derivatives=True,
        derivative_connections=True,
        many_electron_overlaps=True,
    )
    duplicate_snapshot = SimpleNamespace(
        point=SimpleNamespace(
            metadata={CONVERGENCE_METADATA_KEY_V232: convergence.as_dict()}
        ),
        metadata={CONVERGENCE_METADATA_KEY_V232: convergence.as_dict()},
    )
    duplicate_convergence_rejected = _raises_v232(
        lambda: convergence_from_snapshot_v232(duplicate_snapshot),
        (ValueError,),
        "duplicate",
    )

    bundle = build_v231_admission_bundle(
        Path(directory) / "synthetic_live", source_kind="live_ab_initio"
    )
    contract_identity = bundle["provider"].molecular_soc_contract.identity
    method_identity = BackendMethodIdentityV232(
        backend_name=contract_identity.backend_name,
        backend_version=contract_identity.backend_version,
        source_kind=contract_identity.source_kind,
        adapter_name="never-invoked-v232-validator",
        adapter_version="1",
        electronic_method=contract_identity.electronic_method,
        basis=contract_identity.basis,
        active_space=contract_identity.active_space,
        soc_operator=contract_identity.soc_operator,
        scalar_relativistic_method=contract_identity.scalar_relativistic_method,
        derivative_method=contract_identity.derivative_method,
        nac_convention=PYSCF_NAC_CONVENTION_V232,
    ).validate()

    class NeverInvokedValidator:
        parser_name = "never-invoked-parser"
        parser_version = "1"

        def __init__(self):
            self.method_identity = method_identity
            self.calls = 0

        def validate_raw_artifacts_v232(self, **kwargs):
            self.calls += 1
            raise AssertionError("validator executed without prerequisites")

    validator = NeverInvokedValidator()
    policy = BackendAdmissionPolicyV232(
        expected_identity=method_identity,
        trusted_validator_type=NeverInvokedValidator,
        parser_name=validator.parser_name,
        parser_version=validator.parser_version,
    ).validate()
    report = audit_molecular_soc_provider_v232(
        bundle["provider"],
        bundle["center"],
        bundle["dossier_path"],
        requirement="live",
        policy=policy,
        backend_validator=validator,
        execution_challenge="a" * 64,
    )
    inherited_protocol_preserved = bool(
        report.checks["v0231_raw_evidence_protocol"]
    )
    unattested_live_rejected = bool(
        not report.live_backend_admitted
        and not report.checks["runtime_attestation_present"]
    )
    parser_not_invoked = validator.calls == 0

    incomplete_legacy_adapter_rejected = _raises_v232(
        require_pyscf_soc_adapter_v231,
        (ImportError, RuntimeError),
    )
    return {
        "structure_checked_before_capabilities": structure_before_flags,
        "noncallable_engine_methods_rejected": noncallable_rejected,
        "exact_method_identity_enforced": identity_mismatch_rejected,
        "incomplete_capabilities_rejected": incomplete_capabilities_rejected,
        "legacy_convergence_vocabulary_rejected": legacy_convergence_rejected,
        "duplicate_convergence_namespace_rejected": duplicate_convergence_rejected,
        "inherited_raw_evidence_protocol_preserved": inherited_protocol_preserved,
        "unattested_live_source_rejected": unattested_live_rejected,
        "parser_not_invoked_without_prerequisites": parser_not_invoked,
        "incomplete_legacy_pyscf_soc_adapter_rejected": (
            incomplete_legacy_adapter_rejected
        ),
        "unattested_live_report": report.as_dict(),
    }


def run_v0232_release_benchmark():
    inherited = run_v0231_release_benchmark()
    runtime = run_pyscf_runtime_evidence_v232(memory_probe_policy="proc_self")
    thresholds = V232AcceptanceThresholds()
    if len(runtime.checks) != thresholds.expected_runtime_gates:
        raise AssertionError("v0.23.2 runtime evidence must define exactly 28 gates.")

    with tempfile.TemporaryDirectory(prefix="gnd-v232-") as temporary:
        overlap = _overlap_controls_v232(Path(temporary) / "overlap")
        admission = _admission_controls_v232(Path(temporary) / "admission")

    runtime_checks = {
        f"pyscf_runtime::{name}": bool(value)
        for name, value in runtime.checks.items()
    }
    correction_checks = {
        "overlap::physical_contraction_accepted": overlap[
            "physical_contraction_accepted"
        ],
        "overlap::physical_contraction_diagnostics": overlap[
            "physical_contraction_diagnostics"
        ],
        "overlap::spectral_expansion_rejected": overlap[
            "spectral_expansion_rejected"
        ],
        "overlap::reciprocity_corruption_rejected": overlap[
            "reciprocity_corruption_rejected"
        ],
        "overlap::self_identity_corruption_rejected": overlap[
            "self_identity_corruption_rejected"
        ],
        "admission::structure_before_capabilities": admission[
            "structure_checked_before_capabilities"
        ],
        "admission::noncallable_engine_methods_rejected": admission[
            "noncallable_engine_methods_rejected"
        ],
        "admission::exact_method_identity_enforced": admission[
            "exact_method_identity_enforced"
        ],
        "admission::incomplete_capabilities_rejected": admission[
            "incomplete_capabilities_rejected"
        ],
        "admission::legacy_convergence_vocabulary_rejected": admission[
            "legacy_convergence_vocabulary_rejected"
        ],
        "admission::duplicate_convergence_namespace_rejected": admission[
            "duplicate_convergence_namespace_rejected"
        ],
        "admission::inherited_raw_evidence_protocol_preserved": admission[
            "inherited_raw_evidence_protocol_preserved"
        ],
        "admission::unattested_live_source_rejected": admission[
            "unattested_live_source_rejected"
        ],
        "admission::parser_not_invoked_without_prerequisites": admission[
            "parser_not_invoked_without_prerequisites"
        ],
        "admission::incomplete_legacy_soc_adapter_rejected": admission[
            "incomplete_legacy_pyscf_soc_adapter_rejected"
        ],
        "nac_mapping::empirical_mapping_distinct_from_upstream_text": (
            PYSCF_NAC_EMPIRICAL_MAPPING_V232
            != PYSCF_NAC_UPSTREAM_DOCUMENTATION_V232
        ),
        "nac_mapping::internal_indices_pass_through_to_pyscf_tuple": (
            pyscf_state_tuple_for_internal_dij_v232(0, 2) == (0, 2)
        ),
    }
    correction_checks = {
        name: bool(value) for name, value in correction_checks.items()
    }
    new_checks = {**runtime_checks, **correction_checks}
    if len(new_checks) != thresholds.expected_new_gates:
        raise AssertionError("v0.23.2 campaign must define exactly 45 new gates.")

    inherited_checks = {
        f"inherited_v0231::{name}": bool(value)
        for name, value in inherited["acceptance"]["checks"].items()
    }
    if len(inherited_checks) != thresholds.expected_inherited_gates:
        raise AssertionError("v0.23.2 must inherit exactly 123 v0.23.1 gates.")
    checks = {**inherited_checks, **new_checks}
    if len(checks) != thresholds.expected_total_gates:
        raise AssertionError("v0.23.2 campaign must define exactly 168 total gates.")

    return {
        "release": "v0.23.2",
        "theme": (
            "real PySCF spin-free runtime validation, empirical NAC mapping, "
            "finite-manifold overlap correction, and trust-anchored SOC admission"
        ),
        "pyscf_runtime_evidence": runtime.as_dict(),
        "overlap_contract": overlap,
        "runtime_admission_controls": admission,
        "claims": {
            "real_PySCF_spin_free_runtime_validated": True,
            "real_PySCF_SA_CASSCF_gradients_validated": True,
            "real_PySCF_NAC_overlap_consistency_validated": True,
            "finite_manifold_overlap_contract_validated": True,
            "trust_anchored_runtime_admission_validated": True,
            "external_molecular_SOC_snapshot_admitted": False,
            "live_molecular_SOC_backend_admitted": False,
            "ab_initio_SOC_validated": False,
            "live_PySCF_SOC_runtime_validated": False,
            "physical_analytic_SOC_inherited": True,
        },
        "inherited_v0231": inherited,
        "acceptance": {
            "passed": bool(all(checks.values())),
            "checks": checks,
            "inherited_gate_count": len(inherited_checks),
            "runtime_gate_count": len(runtime_checks),
            "new_gate_count": len(new_checks),
            "total_gate_count": len(checks),
            "thresholds": asdict(thresholds),
        },
    }
