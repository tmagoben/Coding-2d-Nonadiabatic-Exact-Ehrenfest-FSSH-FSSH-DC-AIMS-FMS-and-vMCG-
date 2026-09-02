"""Independent accuracy and invariance evidence for v0.24.0 admission."""

from dataclasses import dataclass
import json
from pathlib import Path

import numpy as np

from .openmolcas_rassi_snapshot_v240 import (
    OPENMOLCAS_VALIDATION_NAME_V240,
    ParsedOpenMolcasBundleV240,
    sha256_file_v240,
)


EXTERNAL_VALIDATION_SCHEMA_V240 = "gnd-external-soc-validation-v0.24.0"


def _sha_v240(name, value):
    value = str(value)
    if len(value) != 64 or any(item not in "0123456789abcdef" for item in value):
        raise ValueError(f"{name} must be a lowercase SHA-256 digest.")
    return value


def _matrix_v240(payload, stem, shape):
    real = np.asarray(payload[f"{stem}_real"], dtype=float)
    imag = np.asarray(payload[f"{stem}_imag"], dtype=float)
    if real.shape != shape or imag.shape != shape:
        raise ValueError(f"{stem} has the wrong shape.")
    result = real + 1j * imag
    if not np.all(np.isfinite(result)):
        raise ValueError(f"{stem} contains non-finite values.")
    return result


@dataclass(frozen=True)
class ExternalSOCValidationAuditV240:
    source_kind: str
    reference_error_hartree: float
    basis_changes_hartree: tuple[float, ...]
    method_changes_hartree: tuple[float, ...]
    translation_error_hartree: float
    rotation_error_hartree: float
    checks: dict
    passed: bool

    def as_dict(self):
        return {
            "source_kind": self.source_kind,
            "reference_error_hartree": float(self.reference_error_hartree),
            "basis_changes_hartree": list(self.basis_changes_hartree),
            "method_changes_hartree": list(self.method_changes_hartree),
            "translation_error_hartree": float(self.translation_error_hartree),
            "rotation_error_hartree": float(self.rotation_error_hartree),
            "checks": dict(self.checks),
            "passed": bool(self.passed),
        }


def audit_external_soc_validation_v240(bundle):
    if type(bundle) is not ParsedOpenMolcasBundleV240:
        raise TypeError("external validation requires a parsed OpenMolcas bundle.")
    path = bundle.directory / OPENMOLCAS_VALIDATION_NAME_V240
    if sha256_file_v240(path) != bundle.manifest.validation_sha256:
        raise ValueError("external validation digest differs from the manifest.")
    payload = json.loads(
        path.read_text(encoding="utf-8"),
        parse_constant=lambda value: (_ for _ in ()).throw(
            ValueError(f"non-finite JSON constant {value!r} is forbidden.")
        ),
    )
    expected = {
        "schema",
        "source_kind",
        "protocol_fingerprint",
        "reference_export_sha256",
        "independent_backend_name",
        "independent_backend_version",
        "independent_artifact_sha256",
        "state_order",
        "reference_soc_real",
        "reference_soc_imag",
        "independent_soc_real",
        "independent_soc_imag",
        "reference_tolerance_hartree",
        "basis_labels",
        "basis_soc_real",
        "basis_soc_imag",
        "basis_artifact_sha256",
        "basis_tolerance_hartree",
        "method_labels",
        "method_soc_real",
        "method_soc_imag",
        "method_artifact_sha256",
        "method_tolerance_hartree",
        "frame_base_eigenvalues_hartree",
        "frame_translated_eigenvalues_hartree",
        "frame_rotated_eigenvalues_hartree",
        "frame_artifact_sha256",
        "frame_tolerance_hartree",
        "tracking_minimum_singular_values",
        "tracking_maximum_competing_leakage",
        "tracking_assignment_margins",
        "tracking_artifact_sha256",
    }
    if not isinstance(payload, dict) or set(payload) != expected:
        raise ValueError("external validation fields are incomplete or unknown.")
    if payload["schema"] != EXTERNAL_VALIDATION_SCHEMA_V240:
        raise ValueError("external validation schema mismatch.")
    if payload["source_kind"] != bundle.source_kind:
        raise ValueError("validation and bundle source kinds disagree.")
    if payload["protocol_fingerprint"] != bundle.protocol.fingerprint():
        raise ValueError("validation protocol fingerprint mismatch.")
    reference_record = bundle.manifest.records[0]
    if payload["reference_export_sha256"] != reference_record.export_sha256:
        raise ValueError("validation does not bind the parsed reference export.")
    if tuple(payload["state_order"]) != bundle.protocol.state_order:
        raise ValueError("validation state order mismatch.")
    for name in (
        "independent_backend_name",
        "independent_backend_version",
    ):
        if not str(payload[name]).strip():
            raise ValueError(f"{name} cannot be empty.")
    _sha_v240("independent_artifact_sha256", payload["independent_artifact_sha256"])
    for collection in (
        "basis_artifact_sha256",
        "method_artifact_sha256",
        "frame_artifact_sha256",
        "tracking_artifact_sha256",
    ):
        if not isinstance(payload[collection], list) or not payload[collection]:
            raise ValueError(f"{collection} must be a nonempty artifact digest list.")
        for value in payload[collection]:
            _sha_v240(collection, value)
    nstate = len(bundle.protocol.state_order)
    reference = _matrix_v240(payload, "reference_soc", (nstate, nstate))
    independent = _matrix_v240(payload, "independent_soc", (nstate, nstate))
    parsed_reference = bundle.record_map["reference"].H_soc
    reference_tolerance = float(payload["reference_tolerance_hartree"])
    basis_tolerance = float(payload["basis_tolerance_hartree"])
    method_tolerance = float(payload["method_tolerance_hartree"])
    frame_tolerance = float(payload["frame_tolerance_hartree"])
    if any(
        not np.isfinite(value) or value <= 0.0
        for value in (
            reference_tolerance,
            basis_tolerance,
            method_tolerance,
            frame_tolerance,
        )
    ):
        raise ValueError("validation tolerances must be finite and positive.")
    basis_labels = tuple(str(item).strip() for item in payload["basis_labels"])
    method_labels = tuple(str(item).strip() for item in payload["method_labels"])
    if len(basis_labels) < 3 or len(set(basis_labels)) != len(basis_labels):
        raise ValueError("basis ladder requires at least three unique levels.")
    if len(method_labels) < 3 or len(set(method_labels)) != len(method_labels):
        raise ValueError("method ladder requires at least three unique levels.")
    basis = _matrix_v240(
        payload, "basis_soc", (len(basis_labels), nstate, nstate)
    )
    method = _matrix_v240(
        payload, "method_soc", (len(method_labels), nstate, nstate)
    )
    if len(payload["basis_artifact_sha256"]) != len(basis_labels):
        raise ValueError("basis ladder artifact inventory mismatch.")
    if len(payload["method_artifact_sha256"]) != len(method_labels):
        raise ValueError("method ladder artifact inventory mismatch.")
    basis_changes = tuple(
        float(np.max(np.abs(right - left)))
        for left, right in zip(basis[:-1], basis[1:])
    )
    method_changes = tuple(
        float(np.max(np.abs(right - left)))
        for left, right in zip(method[:-1], method[1:])
    )
    frame_base = np.asarray(payload["frame_base_eigenvalues_hartree"], dtype=float)
    frame_translated = np.asarray(
        payload["frame_translated_eigenvalues_hartree"], dtype=float
    )
    frame_rotated = np.asarray(
        payload["frame_rotated_eigenvalues_hartree"], dtype=float
    )
    if any(item.shape != (nstate,) for item in (frame_base, frame_translated, frame_rotated)):
        raise ValueError("frame spectra must have one value per SOC state.")
    translation_error = float(np.max(np.abs(frame_translated - frame_base)))
    rotation_error = float(np.max(np.abs(frame_rotated - frame_base)))
    minimum_singular = np.asarray(
        payload["tracking_minimum_singular_values"], dtype=float
    )
    leakage = np.asarray(payload["tracking_maximum_competing_leakage"], dtype=float)
    margins = np.asarray(payload["tracking_assignment_margins"], dtype=float)
    expected_tracking = len(bundle.records) - 1
    if any(item.shape != (expected_tracking,) for item in (minimum_singular, leakage, margins)):
        raise ValueError("tracking evidence must cover every displaced record.")
    reference_error = float(np.max(np.abs(reference - independent)))
    checks = {
        "parsed_reference_bound": bool(
            np.max(np.abs(reference - parsed_reference)) <= 1.0e-12
        ),
        "independent_backend": bool(
            payload["source_kind"] == "protocol_fixture"
            or payload["independent_backend_name"].strip().lower() != "openmolcas"
        ),
        "independent_reference_agreement": reference_error <= reference_tolerance,
        "basis_convergence": basis_changes[-1] <= basis_tolerance,
        "method_convergence": method_changes[-1] <= method_tolerance,
        "translation_invariance": translation_error <= frame_tolerance,
        "rotation_invariance": rotation_error <= frame_tolerance,
        "tracking_retention": bool(np.all(minimum_singular >= 0.5)),
        "tracking_leakage": bool(np.all(leakage <= 0.2)),
        "tracking_margin": bool(np.all(margins >= 0.3)),
        "finite_validation_data": bool(
            all(
                np.all(np.isfinite(item))
                for item in (
                    reference,
                    independent,
                    basis,
                    method,
                    frame_base,
                    frame_translated,
                    frame_rotated,
                    minimum_singular,
                    leakage,
                    margins,
                )
            )
        ),
    }
    checks = {name: bool(value) for name, value in checks.items()}
    return ExternalSOCValidationAuditV240(
        source_kind=payload["source_kind"],
        reference_error_hartree=reference_error,
        basis_changes_hartree=basis_changes,
        method_changes_hartree=method_changes,
        translation_error_hartree=translation_error,
        rotation_error_hartree=rotation_error,
        checks=checks,
        passed=bool(all(checks.values())),
    )
