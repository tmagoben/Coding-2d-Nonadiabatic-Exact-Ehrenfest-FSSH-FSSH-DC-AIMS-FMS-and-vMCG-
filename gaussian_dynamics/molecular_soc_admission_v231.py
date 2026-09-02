"""Receipt- and raw-evidence-derived molecular SOC admission for v0.23.1."""

from dataclasses import dataclass
from pathlib import Path
import numpy as np

from .molecular_soc_admission_v230 import audit_molecular_soc_provider_v230
from .molecular_soc_contract_v230 import molecular_soc_contract_from_provider_v230
from .molecular_soc_dossier_v231 import (
    MolecularSOCAdmissionDossierV231,
    load_molecular_soc_dossier_v231,
)


_REQUIREMENTS_V231 = {"protocol", "real", "external", "live"}


@dataclass(frozen=True)
class MolecularSOCAdmissionAuditV231:
    q: np.ndarray
    requirement: str
    source_kind: str
    backend_name: str
    backend_version: str
    replay_dataset_fingerprint: str
    dossier_fingerprint: str
    inherited_v230_report: object
    dossier_protocol_checks: dict
    external_admission_checks: dict
    live_admission_checks: dict
    protocol_passed: bool
    external_snapshot_admitted: bool
    live_backend_admitted: bool
    real_backend_admitted: bool
    passed: bool
    derived_evidence: object
    tracking_report: dict

    def as_dict(self):
        return {
            "q": np.asarray(self.q, dtype=float).tolist(),
            "requirement": self.requirement,
            "source_kind": self.source_kind,
            "backend_name": self.backend_name,
            "backend_version": self.backend_version,
            "replay_dataset_fingerprint": self.replay_dataset_fingerprint,
            "dossier_fingerprint": self.dossier_fingerprint,
            "inherited_v230_report": self.inherited_v230_report.as_dict(),
            "dossier_protocol_checks": {
                name: bool(value)
                for name, value in self.dossier_protocol_checks.items()
            },
            "external_admission_checks": {
                name: bool(value)
                for name, value in self.external_admission_checks.items()
            },
            "live_admission_checks": {
                name: bool(value) for name, value in self.live_admission_checks.items()
            },
            "protocol_passed": bool(self.protocol_passed),
            "external_snapshot_admitted": bool(self.external_snapshot_admitted),
            "live_backend_admitted": bool(self.live_backend_admitted),
            "real_backend_admitted": bool(self.real_backend_admitted),
            "passed": bool(self.passed),
            "derived_evidence": self.derived_evidence.as_dict(),
            "tracking_report": self.tracking_report,
        }


def _resolve_dossier_v231(dossier, dataset, identity, bundle_directory=None):
    if isinstance(dossier, MolecularSOCAdmissionDossierV231):
        if bundle_directory is None:
            raise ValueError(
                "bundle_directory is required when auditing an in-memory dossier."
            )
        return dossier.validate(
            bundle_directory=bundle_directory,
            dataset=dataset,
            identity=identity,
        )
    return load_molecular_soc_dossier_v231(
        Path(dossier), dataset=dataset, identity=identity
    )


def audit_molecular_soc_provider_v231(
    provider,
    q,
    dossier,
    *,
    requirement="real",
    bundle_directory=None,
    backend_validator=None,
):
    """Audit a file-backed admission bundle without trusting summary evidence."""
    if requirement not in _REQUIREMENTS_V231:
        raise ValueError(
            "v0.23.1 admission requirement must be protocol, real, external, or live."
        )
    q = np.asarray(q, dtype=float)
    if q.ndim != 1 or len(q) < 1 or not np.all(np.isfinite(q)):
        raise ValueError("v0.23.1 molecular SOC audit requires a finite coordinate vector.")
    if not hasattr(provider, "dataset"):
        raise TypeError(
            "v0.23.1 admission requires a replay-backed provider with exact raw-data binding."
        )
    dataset = provider.dataset.validate()
    contract = molecular_soc_contract_from_provider_v230(provider)
    identity = contract.identity
    dossier_reference = dossier
    dossier = _resolve_dossier_v231(
        dossier_reference, dataset, identity, bundle_directory=bundle_directory
    )
    if bundle_directory is not None:
        resolved_bundle_directory = Path(bundle_directory)
    elif isinstance(dossier_reference, MolecularSOCAdmissionDossierV231):
        raise ValueError("bundle_directory is required for an in-memory dossier.")
    else:
        dossier_path = Path(dossier_reference)
        resolved_bundle_directory = (
            dossier_path if dossier_path.is_dir() else dossier_path.parent
        )
    inherited = audit_molecular_soc_provider_v230(
        provider, q, require_real_backend=True
    )
    derived = dossier.derived_v230_evidence(dataset)
    tracking = dossier.evidence.tracking.derive(dataset.overlaps)
    receipt_map = {receipt.record_id: receipt for receipt in dossier.receipts}
    trajectory_receipts = [
        receipt_map[record_id] for record_id in dossier.trajectory_record_ids
    ]
    all_receipts_converged = bool(
        dossier.receipts and all(receipt.all_converged for receipt in dossier.receipts)
    )
    trajectory_converged = bool(
        trajectory_receipts
        and all(receipt.all_converged for receipt in trajectory_receipts)
    )
    evidence_matches_contract = bool(
        derived.as_dict() == contract.evidence.as_dict()
    )
    dossier_checks = {
        "inherited_v230_protocol": bool(inherited.protocol_passed),
        "replay_dataset_binding": bool(
            dossier.replay_dataset_fingerprint == dataset.dataset_fingerprint
        ),
        "raw_artifact_integrity": True,
        "trajectory_receipt_coverage": bool(
            len(trajectory_receipts) == len(dataset.q)
        ),
        "trajectory_receipt_convergence": trajectory_converged,
        "all_evidence_calculations_converged": all_receipts_converged,
        "independent_reference_derived": bool(dossier.evidence.reference.passed),
        "basis_convergence_derived": bool(dossier.evidence.basis.passed),
        "method_convergence_derived": bool(dossier.evidence.method.passed),
        "translation_rotation_invariance_derived": bool(
            dossier.evidence.frame.passed
        ),
        "connected_subspace_tracking_derived": bool(tracking["passed"]),
        "derived_evidence_matches_provider_contract": evidence_matches_contract,
    }
    protocol_passed = bool(all(dossier_checks.values()))

    attestation = dossier.runtime_attestation
    attestation_identity = bool(
        attestation is not None
        and attestation.runtime_name.casefold() == identity.backend_name.casefold()
        and attestation.runtime_version == identity.backend_version
        and attestation.environment_sha256 == identity.environment_sha256
    )
    validator_identity = False
    validator_executed = False
    if attestation_identity and backend_validator is not None:
        validator_identity = bool(
            str(getattr(backend_validator, "adapter_name", ""))
            == attestation.adapter_name
            and str(getattr(backend_validator, "adapter_version", ""))
            == attestation.adapter_version
            and callable(getattr(backend_validator, "validate_raw_artifacts", None))
        )
        if validator_identity:
            result = backend_validator.validate_raw_artifacts(
                dossier=dossier,
                bundle_directory=resolved_bundle_directory,
                dataset=dataset,
            )
            validator_executed = result is True
    external_checks = {
        "v0231_protocol": protocol_passed,
        "inherited_v230_real_admission": bool(inherited.real_backend_admitted),
        "external_snapshot_source": identity.source_kind
        == "external_ab_initio_snapshot",
        "backend_attestation_present": attestation is not None,
        "backend_attestation_identity": attestation_identity,
        "backend_validator_identity": validator_identity,
        "executable_artifact_validation": validator_executed,
        "backend_artifact_parser_validated": bool(
            attestation is not None and attestation.artifact_parser_validated
        ),
        "method_specific_SOC_validated": bool(
            attestation is not None and attestation.external_ready
        ),
    }
    live_checks = {
        "v0231_protocol": protocol_passed,
        "inherited_v230_real_admission": bool(inherited.real_backend_admitted),
        "live_ab_initio_source": identity.source_kind == "live_ab_initio",
        "backend_attestation_present": attestation is not None,
        "backend_attestation_identity": attestation_identity,
        "backend_validator_identity": validator_identity,
        "executable_artifact_validation": validator_executed,
        "backend_artifact_parser_validated": bool(
            attestation is not None and attestation.artifact_parser_validated
        ),
        "fresh_runtime_execution_validated": bool(
            attestation is not None and attestation.live_ready
        ),
    }
    external_admitted = bool(all(external_checks.values()))
    live_admitted = bool(all(live_checks.values()))
    real_admitted = bool(external_admitted or live_admitted)
    passed_by_requirement = {
        "protocol": protocol_passed,
        "real": real_admitted,
        "external": external_admitted,
        "live": live_admitted,
    }
    return MolecularSOCAdmissionAuditV231(
        q=q.copy(),
        requirement=requirement,
        source_kind=identity.source_kind,
        backend_name=identity.backend_name,
        backend_version=identity.backend_version,
        replay_dataset_fingerprint=dataset.dataset_fingerprint,
        dossier_fingerprint=dossier.fingerprint(),
        inherited_v230_report=inherited,
        dossier_protocol_checks={
            name: bool(value) for name, value in dossier_checks.items()
        },
        external_admission_checks={
            name: bool(value) for name, value in external_checks.items()
        },
        live_admission_checks={
            name: bool(value) for name, value in live_checks.items()
        },
        protocol_passed=protocol_passed,
        external_snapshot_admitted=external_admitted,
        live_backend_admitted=live_admitted,
        real_backend_admitted=real_admitted,
        passed=passed_by_requirement[requirement],
        derived_evidence=derived,
        tracking_report=tracking,
    )


def require_molecular_soc_protocol_v231(*args, **kwargs):
    kwargs["requirement"] = "protocol"
    report = audit_molecular_soc_provider_v231(*args, **kwargs)
    if not report.protocol_passed:
        failed = ", ".join(
            name for name, value in report.dossier_protocol_checks.items() if not value
        )
        raise ValueError(f"v0.23.1 molecular SOC protocol failed: {failed}.")
    return report


def require_external_molecular_soc_snapshot_v231(*args, **kwargs):
    kwargs["requirement"] = "external"
    report = audit_molecular_soc_provider_v231(*args, **kwargs)
    if not report.external_snapshot_admitted:
        failed = ", ".join(
            name for name, value in report.external_admission_checks.items() if not value
        )
        raise ValueError(f"external molecular SOC snapshot admission failed: {failed}.")
    return report


def require_live_molecular_soc_backend_v231(*args, **kwargs):
    kwargs["requirement"] = "live"
    report = audit_molecular_soc_provider_v231(*args, **kwargs)
    if not report.live_backend_admitted:
        failed = ", ".join(
            name for name, value in report.live_admission_checks.items() if not value
        )
        raise ValueError(f"live molecular SOC backend admission failed: {failed}.")
    return report
