"""Trust-anchored external and live molecular SOC admission for v0.23.2."""

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .molecular_soc_admission_v231 import audit_molecular_soc_provider_v231
from .molecular_soc_contract_v230 import molecular_soc_contract_from_provider_v230
from .molecular_soc_dossier_v231 import (
    MolecularSOCAdmissionDossierV231,
    load_molecular_soc_dossier_v231,
)
from .molecular_soc_runtime_v232 import (
    BackendAdmissionPolicyV232,
    BackendArtifactValidationProofV232,
    BackendMethodIdentityV232,
    load_runtime_probe_v232,
)


_REQUIREMENTS_V232 = {"external", "live"}


@dataclass(frozen=True)
class MolecularSOCRuntimeAdmissionAuditV232:
    requirement: str
    source_kind: str
    backend_name: str
    backend_version: str
    dossier_fingerprint: str
    execution_challenge: str
    checks: dict
    validator_error: str | None
    external_snapshot_admitted: bool
    live_backend_admitted: bool
    passed: bool

    def as_dict(self):
        return {
            "requirement": self.requirement,
            "source_kind": self.source_kind,
            "backend_name": self.backend_name,
            "backend_version": self.backend_version,
            "dossier_fingerprint": self.dossier_fingerprint,
            "execution_challenge": self.execution_challenge,
            "checks": {name: bool(value) for name, value in self.checks.items()},
            "validator_error": self.validator_error,
            "external_snapshot_admitted": bool(self.external_snapshot_admitted),
            "live_backend_admitted": bool(self.live_backend_admitted),
            "passed": bool(self.passed),
        }


def _resolve_dossier_v232(dossier, dataset, identity, bundle_directory):
    if isinstance(dossier, MolecularSOCAdmissionDossierV231):
        if bundle_directory is None:
            raise ValueError("bundle_directory is required for an in-memory dossier.")
        return (
            dossier.validate(
                bundle_directory=bundle_directory,
                dataset=dataset,
                identity=identity,
            ),
            Path(bundle_directory),
        )
    path = Path(dossier)
    loaded = load_molecular_soc_dossier_v231(
        path, dataset=dataset, identity=identity
    )
    return loaded, path if path.is_dir() else path.parent


def _proof_checks_v232(
    proof,
    *,
    policy,
    identity,
    dossier,
    artifact_map,
    receipt_map,
    challenge,
    runtime_probe,
):
    if type(proof) is not BackendArtifactValidationProofV232:
        return {"typed_parser_execution_proof": False}
    proof.validate()
    parsed_outputs = dict(proof.parsed_output_artifacts)
    expected_outputs = {
        receipt.output_artifact: artifact_map[receipt.output_artifact].sha256
        for receipt in receipt_map.values()
    }
    execution_receipts = {item.record_id: item for item in proof.receipt_evidence}
    exact_receipts = bool(set(execution_receipts) == set(receipt_map))
    if exact_receipts:
        exact_receipts = all(
            execution_receipts[record_id].output_artifact
            == receipt_map[record_id].output_artifact
            for record_id in receipt_map
        )
    convergence_complete = bool(
        exact_receipts
        and all(item.convergence.complete for item in execution_receipts.values())
    )
    attestation = dossier.runtime_attestation
    probe_record = artifact_map[attestation.runtime_probe_artifact]
    return {
        "typed_parser_execution_proof": True,
        "proof_method_identity": proof.method_identity == policy.expected_identity,
        "proof_matches_contract_identity": proof.method_identity.matches_backend_identity(
            identity
        ),
        "trusted_parser_identity": bool(
            proof.parser_name == policy.parser_name
            and proof.parser_version == policy.parser_version
        ),
        "fresh_execution_challenge_bound": proof.execution_challenge == challenge,
        "replay_dataset_bound": proof.replay_dataset_fingerprint
        == dossier.replay_dataset_fingerprint,
        "dossier_fingerprint_bound": proof.dossier_fingerprint
        == dossier.fingerprint(),
        "environment_bound": proof.environment_sha256
        == identity.environment_sha256,
        "calculation_template_bound": proof.calculation_input_sha256
        == identity.calculation_input_sha256,
        "runtime_probe_name_bound": proof.runtime_probe_artifact
        == attestation.runtime_probe_artifact,
        "runtime_probe_hash_bound": proof.runtime_probe_sha256
        == probe_record.sha256,
        "exact_parsed_output_inventory": parsed_outputs == expected_outputs,
        "exact_receipt_execution_inventory": exact_receipts,
        "component_resolved_convergence": convergence_complete,
        "parser_actually_invoked": bool(proof.parser_executed),
        "fresh_backend_execution_observed": bool(proof.fresh_execution_observed),
        "probe_method_identity": runtime_probe.method_identity
        == policy.expected_identity,
        "probe_environment_bound": runtime_probe.environment_sha256
        == identity.environment_sha256,
        "probe_calculation_template_bound": runtime_probe.calculation_input_sha256
        == identity.calculation_input_sha256,
        "probe_replay_dataset_bound": runtime_probe.replay_dataset_fingerprint
        == dossier.replay_dataset_fingerprint,
        "probe_runtime_imported": bool(runtime_probe.runtime_imported),
    }


def audit_molecular_soc_provider_v232(
    provider,
    q,
    dossier,
    *,
    requirement,
    policy,
    backend_validator,
    execution_challenge,
    bundle_directory=None,
):
    """Require independent trust policy plus typed parser/execution evidence.

    ``policy`` is intentionally supplied by the caller.  Reading the trusted adapter
    type or parser identity from the dossier/engine would allow a synthetic engine to
    certify itself.
    """
    if requirement not in _REQUIREMENTS_V232:
        raise ValueError("v0.23.2 requirement must be external or live.")
    if not isinstance(policy, BackendAdmissionPolicyV232):
        raise TypeError("policy must be BackendAdmissionPolicyV232.")
    policy.validate()
    if (
        not isinstance(execution_challenge, str)
        or len(execution_challenge) != 64
        or any(character not in "0123456789abcdef" for character in execution_challenge)
    ):
        raise ValueError("execution_challenge must be a lowercase SHA-256 nonce.")
    q = np.asarray(q, dtype=float)
    if q.ndim != 1 or len(q) < 1 or not np.all(np.isfinite(q)):
        raise ValueError("v0.23.2 audit requires a finite coordinate vector.")
    if not hasattr(provider, "dataset"):
        raise TypeError("v0.23.2 admission requires an exactly replay-bound provider.")

    dataset = provider.dataset.validate()
    contract = molecular_soc_contract_from_provider_v230(provider)
    identity = contract.identity
    resolved_dossier, resolved_directory = _resolve_dossier_v232(
        dossier, dataset, identity, bundle_directory
    )
    inherited = audit_molecular_soc_provider_v231(
        provider,
        q,
        resolved_dossier,
        requirement="protocol",
        bundle_directory=resolved_directory,
    )
    artifact_map, receipt_map = resolved_dossier._maps(resolved_directory)
    attestation = resolved_dossier.runtime_attestation

    checks = {
        "v0231_raw_evidence_protocol": bool(inherited.protocol_passed),
        "policy_matches_contract_identity": policy.expected_identity.matches_backend_identity(
            identity
        ),
        "trusted_validator_exact_type": type(backend_validator)
        is policy.trusted_validator_type,
        "validator_method_identity": bool(
            type(getattr(backend_validator, "method_identity", None))
            is BackendMethodIdentityV232
            and backend_validator.method_identity == policy.expected_identity
        ),
        "validator_parser_callable": callable(
            getattr(backend_validator, "validate_raw_artifacts_v232", None)
        ),
        "runtime_attestation_present": attestation is not None,
        "runtime_attestation_exact_identity": bool(
            attestation is not None
            and attestation.runtime_name == policy.expected_identity.backend_name
            and attestation.runtime_version
            == policy.expected_identity.backend_version
            and attestation.adapter_name == policy.expected_identity.adapter_name
            and attestation.adapter_version
            == policy.expected_identity.adapter_version
            and attestation.environment_sha256 == identity.environment_sha256
        ),
        "runtime_attestation_capabilities": bool(
            attestation is not None and attestation.external_ready
        ),
        "runtime_attestation_live_execution": bool(
            attestation is not None and attestation.live_ready
        ),
    }

    validator_error = None
    proof = None
    runtime_probe = None
    prerequisite_names = set(checks)
    if requirement == "external":
        prerequisite_names.discard("runtime_attestation_live_execution")
    prerequisites = bool(all(checks[name] for name in prerequisite_names))
    if prerequisites:
        try:
            probe_record = artifact_map[attestation.runtime_probe_artifact]
            runtime_probe = load_runtime_probe_v232(
                resolved_directory / probe_record.relative_path
            )
            proof = backend_validator.validate_raw_artifacts_v232(
                dossier=resolved_dossier,
                bundle_directory=resolved_directory,
                dataset=dataset,
                execution_challenge=execution_challenge,
            )
            checks.update(
                _proof_checks_v232(
                    proof,
                    policy=policy,
                    identity=identity,
                    dossier=resolved_dossier,
                    artifact_map=artifact_map,
                    receipt_map=receipt_map,
                    challenge=execution_challenge,
                    runtime_probe=runtime_probe,
                )
            )
        except Exception as exc:  # A backend parser is an untrusted boundary.
            validator_error = f"{type(exc).__name__}: {exc}"
            checks["typed_parser_execution_proof"] = False
    else:
        checks["typed_parser_execution_proof"] = False

    negative_control = identity.extra.get("synthetic_relabel_negative_control")
    molecular_claim = identity.extra.get("molecular_accuracy_claim")
    checks["explicit_negative_control_absent"] = bool(
        negative_control is None
        or isinstance(negative_control, (bool, np.bool_))
        and not bool(negative_control)
    )
    checks["explicit_false_molecular_claim_absent"] = bool(
        molecular_claim is None
        or isinstance(molecular_claim, (bool, np.bool_))
        and bool(molecular_claim)
    )

    source_external = identity.source_kind == "external_ab_initio_snapshot"
    source_live = identity.source_kind == "live_ab_initio"
    checks["external_snapshot_source"] = source_external
    checks["live_ab_initio_source"] = source_live

    shared_names = {
        name
        for name in checks
        if name not in {"external_snapshot_source", "live_ab_initio_source"}
    }
    external_names = shared_names - {
        "fresh_backend_execution_observed",
        "runtime_attestation_live_execution",
    }
    external_admitted = bool(
        source_external and all(checks[name] for name in external_names)
    )
    live_admitted = bool(
        source_live and all(checks[name] for name in shared_names)
    )
    return MolecularSOCRuntimeAdmissionAuditV232(
        requirement=requirement,
        source_kind=identity.source_kind,
        backend_name=identity.backend_name,
        backend_version=identity.backend_version,
        dossier_fingerprint=resolved_dossier.fingerprint(),
        execution_challenge=execution_challenge,
        checks=checks,
        validator_error=validator_error,
        external_snapshot_admitted=external_admitted,
        live_backend_admitted=live_admitted,
        passed=external_admitted if requirement == "external" else live_admitted,
    )


def require_external_molecular_soc_snapshot_v232(*args, **kwargs):
    kwargs["requirement"] = "external"
    report = audit_molecular_soc_provider_v232(*args, **kwargs)
    if not report.external_snapshot_admitted:
        failed = ", ".join(name for name, value in report.checks.items() if not value)
        raise ValueError("v0.23.2 external SOC admission failed: " + failed)
    return report


def require_live_molecular_soc_backend_v232(*args, **kwargs):
    kwargs["requirement"] = "live"
    report = audit_molecular_soc_provider_v232(*args, **kwargs)
    if not report.live_backend_admitted:
        failed = ", ".join(name for name, value in report.checks.items() if not value)
        raise ValueError("v0.23.2 live SOC admission failed: " + failed)
    return report
