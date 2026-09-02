"""Trust-anchored admission of an external OpenMolcas SOC snapshot."""

from dataclasses import dataclass

import numpy as np

from .external_soc_validation_v240 import audit_external_soc_validation_v240
from .molecular_soc_convention_v233 import (
    MolecularSOCMatrixConventionV233,
    require_exact_molecular_soc_convention_v233,
)
from .openmolcas_rassi_protocol_v240 import OpenMolcasRASSIProtocolV240
from .openmolcas_rassi_snapshot_v240 import (
    NATIVE_OPENMOLCAS_NUMERIC_CROSSCHECK_V240,
    OpenMolcasRASSISnapshotParserV240,
    ParsedOpenMolcasBundleV240,
)
from .soc_admission_v221 import audit_soc_symmetry_contract_v221
from .soc_derivative_evidence_v240 import audit_external_soc_derivatives_v240


def _sha_v240(name, value):
    value = str(value)
    if len(value) != 64 or any(item not in "0123456789abcdef" for item in value):
        raise ValueError(f"{name} must be a lowercase SHA-256 digest.")
    return value


@dataclass(frozen=True)
class ExternalSOCAdmissionPolicyV240:
    """Out-of-band trust policy; never derive this from the submitted bundle."""

    expected_protocol: OpenMolcasRASSIProtocolV240
    expected_soc_convention: MolecularSOCMatrixConventionV233
    trusted_parser_type: type
    expected_manifest_sha256: str
    expected_environment_sha256: str
    expected_exporter_name: str
    expected_exporter_version: str

    def validate(self):
        self.expected_protocol.validate()
        self.expected_soc_convention.validate()
        if self.expected_soc_convention != self.expected_protocol.soc_convention():
            raise ValueError("trusted protocol and SOC convention disagree.")
        if self.trusted_parser_type is not OpenMolcasRASSISnapshotParserV240:
            raise TypeError("v0.24.0 trusts only the exact OpenMolcas parser type.")
        _sha_v240("expected_manifest_sha256", self.expected_manifest_sha256)
        _sha_v240("expected_environment_sha256", self.expected_environment_sha256)
        if not self.expected_exporter_name.strip() or not self.expected_exporter_version.strip():
            raise ValueError("trusted exporter identity cannot be empty.")
        return self


@dataclass(frozen=True)
class ExternalSOCAdmissionAuditV240:
    source_kind: str
    protocol_fingerprint: str
    bundle_fingerprint: str | None
    parser_error: str | None
    checks: dict
    protocol_passed: bool
    external_snapshot_admitted: bool
    live_backend_admitted: bool

    @property
    def passed(self):
        return self.external_snapshot_admitted

    def as_dict(self):
        return {
            "source_kind": self.source_kind,
            "protocol_fingerprint": self.protocol_fingerprint,
            "bundle_fingerprint": self.bundle_fingerprint,
            "parser_error": self.parser_error,
            "checks": dict(self.checks),
            "protocol_passed": bool(self.protocol_passed),
            "external_snapshot_admitted": bool(self.external_snapshot_admitted),
            "live_backend_admitted": False,
            "passed": bool(self.passed),
        }


def _time_reversal_residual_v240(matrix, time_reversal_matrix):
    matrix = np.asarray(matrix, dtype=complex)
    J = np.asarray(time_reversal_matrix, dtype=complex)
    return float(np.linalg.norm(matrix - J @ matrix.conj() @ J.conj().T, ord="fro"))


def audit_external_soc_snapshot_v240(directory, *, policy, parser):
    if type(policy) is not ExternalSOCAdmissionPolicyV240:
        raise TypeError("policy must be ExternalSOCAdmissionPolicyV240.")
    policy.validate()
    checks = {
        "trusted_parser_exact_type": type(parser) is policy.trusted_parser_type,
        "trusted_parser_identity": bool(
            getattr(parser, "parser_name", None) == "gnd-openmolcas-rassi-snapshot"
            and getattr(parser, "parser_version", None) == "0.24.0"
        ),
    }
    parser_error = None
    bundle = None
    derivative = None
    validation = None
    if all(checks.values()):
        try:
            bundle = parser.parse_bundle_v240(directory)
            if type(bundle) is not ParsedOpenMolcasBundleV240:
                raise TypeError("trusted parser returned an unexpected proof type.")
            derivative = audit_external_soc_derivatives_v240(bundle)
            validation = audit_external_soc_validation_v240(bundle)
        except Exception as exc:
            parser_error = f"{type(exc).__name__}: {exc}"
    checks["typed_parser_execution"] = bool(
        bundle is not None and bundle.parser_executed
    )
    if bundle is None:
        checks = {name: bool(value) for name, value in checks.items()}
        return ExternalSOCAdmissionAuditV240(
            source_kind="unparsed",
            protocol_fingerprint=policy.expected_protocol.fingerprint(),
            bundle_fingerprint=None,
            parser_error=parser_error,
            checks=checks,
            protocol_passed=False,
            external_snapshot_admitted=False,
            live_backend_admitted=False,
        )

    protocol = bundle.protocol
    reference = bundle.record_map["reference"]
    symmetry = protocol.symmetry_contract()
    symmetry_audit = audit_soc_symmetry_contract_v221(
        protocol.model_space(), symmetry
    )
    convention_exact = True
    try:
        require_exact_molecular_soc_convention_v233(
            protocol.soc_convention(), policy.expected_soc_convention
        )
    except (TypeError, ValueError):
        convention_exact = False
    tr_residual = max(
        _time_reversal_residual_v240(reference.H_spin_free, symmetry.time_reversal_matrix),
        _time_reversal_residual_v240(reference.H_soc, symmetry.time_reversal_matrix),
        _time_reversal_residual_v240(reference.H_total, symmetry.time_reversal_matrix),
    )
    hermiticity_residual = max(
        np.linalg.norm(reference.H_spin_free - reference.H_spin_free.conj().T, ord="fro"),
        np.linalg.norm(reference.H_soc - reference.H_soc.conj().T, ord="fro"),
    )
    checks.update(
        {
            "exact_protocol_identity": protocol == policy.expected_protocol,
            "exact_protocol_fingerprint": (
                protocol.fingerprint() == policy.expected_protocol.fingerprint()
            ),
            "exact_soc_convention": convention_exact,
            "exact_manifest_digest": (
                bundle.manifest_sha256 == policy.expected_manifest_sha256
            ),
            "exact_environment_digest": (
                bundle.manifest.environment_sha256
                == policy.expected_environment_sha256
            ),
            "exact_exporter_identity": bool(
                bundle.manifest.exporter_name == policy.expected_exporter_name
                and bundle.manifest.exporter_version
                == policy.expected_exporter_version
            ),
            "exact_artifact_inventory": bundle.exact_artifact_inventory,
            "complete_even_electron_model_space": symmetry_audit.passed,
            "reference_component_hermiticity": hermiticity_residual <= 1.0e-10,
            "reference_time_reversal_symmetry": tr_residual <= 1.0e-10,
            "nonzero_soc_signal": np.linalg.norm(reference.H_soc, ord="fro") > 0.0,
            "displaced_derivative_evidence": bool(
                derivative is not None and derivative.passed
            ),
            "independent_accuracy_evidence": bool(
                validation is not None and validation.passed
            ),
            "external_snapshot_source": (
                bundle.source_kind == "external_ab_initio_snapshot"
            ),
            "native_openmolcas_execution": bundle.native_openmolcas_execution,
            "native_numeric_crosscheck": NATIVE_OPENMOLCAS_NUMERIC_CROSSCHECK_V240,
            "live_execution_source": False,
        }
    )
    checks = {name: bool(value) for name, value in checks.items()}
    protocol_only_exclusions = {
        "external_snapshot_source",
        "native_openmolcas_execution",
        "native_numeric_crosscheck",
        "live_execution_source",
    }
    protocol_passed = bool(
        all(value for name, value in checks.items() if name not in protocol_only_exclusions)
    )
    external_checks = {
        name: value for name, value in checks.items() if name != "live_execution_source"
    }
    external_admitted = bool(all(external_checks.values()))
    return ExternalSOCAdmissionAuditV240(
        source_kind=bundle.source_kind,
        protocol_fingerprint=protocol.fingerprint(),
        bundle_fingerprint=bundle.fingerprint,
        parser_error=parser_error,
        checks=checks,
        protocol_passed=protocol_passed,
        external_snapshot_admitted=external_admitted,
        live_backend_admitted=False,
    )


def require_external_soc_snapshot_v240(*args, **kwargs):
    audit = audit_external_soc_snapshot_v240(*args, **kwargs)
    if not audit.external_snapshot_admitted:
        failed = ", ".join(name for name, value in audit.checks.items() if not value)
        raise ValueError("external OpenMolcas SOC snapshot not admitted: " + failed)
    return audit
