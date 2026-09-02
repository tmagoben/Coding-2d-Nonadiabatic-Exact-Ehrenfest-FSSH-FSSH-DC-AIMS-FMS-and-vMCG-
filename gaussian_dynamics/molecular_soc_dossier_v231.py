"""Traceable calculation receipts and admission dossiers for v0.23.1."""

from dataclasses import dataclass
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import numpy as np

from .molecular_soc_contract_v230 import MolecularSOCBackendIdentityV230
from .molecular_soc_evidence_v231 import DerivedEvidenceBundleV231


DOSSIER_NAME_V231 = "molecular_soc_admission_dossier_v231.json"
_SOURCE_KINDS_V231 = {
    "validation_fixture",
    "external_ab_initio_snapshot",
    "live_ab_initio",
}
_RECEIPT_ROLES_V231 = {
    "trajectory",
    "basis",
    "method",
    "frame_base",
    "frame_translation",
    "frame_rotation",
}
_ARTIFACT_ROLES_V231 = {
    "calculation_template",
    "environment_lock",
    "calculation_input",
    "calculation_output",
    "independent_reference",
    "runtime_probe",
}


def _required_text_v231(name, value):
    text = str(value).strip()
    if not text:
        raise ValueError(f"{name} cannot be empty.")
    return text


def _native_bool_v231(name, value):
    if not isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{name} must be Boolean.")
    return bool(value)


def _sha256_digest_v231(name, value):
    value = str(value)
    if len(value) != 64 or any(character not in "0123456789abcdef" for character in value):
        raise ValueError(f"{name} must be a lowercase SHA-256 digest.")
    return value


def _sha256_file_v231(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _canonical_v231(value):
    if isinstance(value, np.generic):
        return _canonical_v231(value.item())
    if isinstance(value, complex):
        if not np.isfinite(value.real) or not np.isfinite(value.imag):
            raise ValueError("dossier cannot contain non-finite complex values.")
        return {"real": float(value.real), "imag": float(value.imag)}
    if isinstance(value, np.ndarray):
        return _canonical_v231(value.tolist())
    if isinstance(value, dict):
        if any(not isinstance(key, str) for key in value):
            raise TypeError("dossier dictionary keys must be strings.")
        return {
            key: _canonical_v231(item) for key, item in sorted(value.items())
        }
    if isinstance(value, (tuple, list)):
        return [_canonical_v231(item) for item in value]
    if isinstance(value, float) and not np.isfinite(value):
        raise ValueError("dossier cannot contain non-finite values.")
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    raise TypeError(f"unsupported dossier value {type(value).__name__}.")


def _canonical_json_bytes_v231(value):
    return json.dumps(
        _canonical_v231(value), sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")


@dataclass(frozen=True)
class RawArtifactRecordV231:
    name: str
    relative_path: str
    role: str
    sha256: str
    size_bytes: int

    def validate(self, bundle_directory=None):
        _required_text_v231("artifact name", self.name)
        relative = PurePosixPath(str(self.relative_path))
        if relative.is_absolute() or not relative.parts or ".." in relative.parts:
            raise ValueError("artifact relative_path must stay inside the dossier bundle.")
        if self.role not in _ARTIFACT_ROLES_V231:
            raise ValueError(f"unsupported raw artifact role {self.role!r}.")
        _sha256_digest_v231("artifact sha256", self.sha256)
        if int(self.size_bytes) != self.size_bytes or int(self.size_bytes) < 0:
            raise ValueError("artifact size_bytes must be a nonnegative integer.")
        if bundle_directory is not None:
            root = Path(bundle_directory).resolve()
            path = (root / Path(*relative.parts)).resolve()
            if not path.is_relative_to(root) or not path.is_file():
                raise ValueError(f"raw artifact is missing: {self.relative_path}.")
            if path.stat().st_size != int(self.size_bytes):
                raise ValueError(f"raw artifact size mismatch: {self.name}.")
            if _sha256_file_v231(path) != self.sha256:
                raise ValueError(f"raw artifact integrity mismatch: {self.name}.")
        return self

    def as_dict(self):
        self.validate()
        return {
            "name": self.name,
            "relative_path": self.relative_path,
            "role": self.role,
            "sha256": self.sha256,
            "size_bytes": int(self.size_bytes),
        }

    @classmethod
    def from_dict(cls, payload):
        return cls(
            name=payload["name"],
            relative_path=payload["relative_path"],
            role=payload["role"],
            sha256=payload["sha256"],
            size_bytes=payload["size_bytes"],
        ).validate()


def write_raw_json_artifact_v231(
    bundle_directory,
    *,
    name,
    relative_path,
    role,
    payload,
    overwrite=False,
):
    """Write one canonical raw artifact and return its immutable record."""
    bundle_directory = Path(bundle_directory)
    record_path = bundle_directory / relative_path
    record_path.parent.mkdir(parents=True, exist_ok=True)
    if record_path.exists() and not overwrite:
        raise FileExistsError(f"raw artifact already exists: {relative_path}.")
    temporary = record_path.parent / f".{record_path.name}.tmp"
    temporary.write_bytes(_canonical_json_bytes_v231(payload))
    os.replace(temporary, record_path)
    return RawArtifactRecordV231(
        name=name,
        relative_path=PurePosixPath(relative_path).as_posix(),
        role=role,
        sha256=_sha256_file_v231(record_path),
        size_bytes=record_path.stat().st_size,
    ).validate(bundle_directory)


@dataclass(frozen=True)
class CalculationReceiptV231:
    record_id: str
    role: str
    q_bohr: tuple[float, ...]
    backend_name: str
    backend_version: str
    source_kind: str
    electronic_method: str
    basis: str
    soc_operator: str
    derivative_method: str
    input_artifact: str
    output_artifact: str
    scf_converged: bool
    correlated_converged: bool
    soc_converged: bool
    derivatives_converged: bool
    overlaps_converged: bool

    def validate(self):
        for name in (
            "record_id",
            "backend_name",
            "backend_version",
            "electronic_method",
            "basis",
            "soc_operator",
            "derivative_method",
            "input_artifact",
            "output_artifact",
        ):
            _required_text_v231(name, getattr(self, name))
        if self.role not in _RECEIPT_ROLES_V231:
            raise ValueError(f"unsupported calculation receipt role {self.role!r}.")
        if self.source_kind not in _SOURCE_KINDS_V231:
            raise ValueError(f"unsupported calculation source kind {self.source_kind!r}.")
        q = np.asarray(self.q_bohr, dtype=float)
        if q.ndim != 1 or len(q) < 1 or not np.all(np.isfinite(q)):
            raise ValueError("calculation receipt q_bohr must be a finite vector.")
        for name in (
            "scf_converged",
            "correlated_converged",
            "soc_converged",
            "derivatives_converged",
            "overlaps_converged",
        ):
            _native_bool_v231(name, getattr(self, name))
        if self.input_artifact == self.output_artifact:
            raise ValueError("calculation input and output artifacts must differ.")
        return self

    @property
    def all_converged(self):
        self.validate()
        return bool(
            self.scf_converged
            and self.correlated_converged
            and self.soc_converged
            and self.derivatives_converged
            and self.overlaps_converged
        )

    def as_dict(self):
        self.validate()
        return {
            "record_id": self.record_id,
            "role": self.role,
            "q_bohr": [float(value) for value in self.q_bohr],
            "backend_name": self.backend_name,
            "backend_version": self.backend_version,
            "source_kind": self.source_kind,
            "electronic_method": self.electronic_method,
            "basis": self.basis,
            "soc_operator": self.soc_operator,
            "derivative_method": self.derivative_method,
            "input_artifact": self.input_artifact,
            "output_artifact": self.output_artifact,
            "scf_converged": bool(self.scf_converged),
            "correlated_converged": bool(self.correlated_converged),
            "soc_converged": bool(self.soc_converged),
            "derivatives_converged": bool(self.derivatives_converged),
            "overlaps_converged": bool(self.overlaps_converged),
            "all_converged": self.all_converged,
        }

    @classmethod
    def from_dict(cls, payload):
        return cls(
            record_id=payload["record_id"],
            role=payload["role"],
            q_bohr=tuple(payload["q_bohr"]),
            backend_name=payload["backend_name"],
            backend_version=payload["backend_version"],
            source_kind=payload["source_kind"],
            electronic_method=payload["electronic_method"],
            basis=payload["basis"],
            soc_operator=payload["soc_operator"],
            derivative_method=payload["derivative_method"],
            input_artifact=payload["input_artifact"],
            output_artifact=payload["output_artifact"],
            scf_converged=payload["scf_converged"],
            correlated_converged=payload["correlated_converged"],
            soc_converged=payload["soc_converged"],
            derivatives_converged=payload["derivatives_converged"],
            overlaps_converged=payload["overlaps_converged"],
        ).validate()


@dataclass(frozen=True)
class BackendRuntimeAttestationV231:
    runtime_name: str
    runtime_version: str
    adapter_name: str
    adapter_version: str
    environment_sha256: str
    runtime_probe_artifact: str
    runtime_imported: bool
    method_specific_soc_implemented: bool
    soc_derivatives_implemented: bool
    wavefunction_overlaps_implemented: bool
    artifact_parser_validated: bool
    fresh_execution_observed: bool

    def validate(self):
        for name in (
            "runtime_name",
            "runtime_version",
            "adapter_name",
            "adapter_version",
            "runtime_probe_artifact",
        ):
            _required_text_v231(name, getattr(self, name))
        _sha256_digest_v231("runtime environment_sha256", self.environment_sha256)
        for name in (
            "runtime_imported",
            "method_specific_soc_implemented",
            "soc_derivatives_implemented",
            "wavefunction_overlaps_implemented",
            "artifact_parser_validated",
            "fresh_execution_observed",
        ):
            _native_bool_v231(name, getattr(self, name))
        return self

    @property
    def external_ready(self):
        self.validate()
        return bool(
            self.method_specific_soc_implemented
            and self.soc_derivatives_implemented
            and self.wavefunction_overlaps_implemented
            and self.artifact_parser_validated
        )

    @property
    def live_ready(self):
        return bool(
            self.external_ready
            and self.runtime_imported
            and self.fresh_execution_observed
        )

    def as_dict(self):
        self.validate()
        return {
            "runtime_name": self.runtime_name,
            "runtime_version": self.runtime_version,
            "adapter_name": self.adapter_name,
            "adapter_version": self.adapter_version,
            "environment_sha256": self.environment_sha256,
            "runtime_probe_artifact": self.runtime_probe_artifact,
            "runtime_imported": bool(self.runtime_imported),
            "method_specific_soc_implemented": bool(
                self.method_specific_soc_implemented
            ),
            "soc_derivatives_implemented": bool(self.soc_derivatives_implemented),
            "wavefunction_overlaps_implemented": bool(
                self.wavefunction_overlaps_implemented
            ),
            "artifact_parser_validated": bool(self.artifact_parser_validated),
            "fresh_execution_observed": bool(self.fresh_execution_observed),
            "external_ready": self.external_ready,
            "live_ready": self.live_ready,
        }

    @classmethod
    def from_dict(cls, payload):
        return cls(
            runtime_name=payload["runtime_name"],
            runtime_version=payload["runtime_version"],
            adapter_name=payload["adapter_name"],
            adapter_version=payload["adapter_version"],
            environment_sha256=payload["environment_sha256"],
            runtime_probe_artifact=payload["runtime_probe_artifact"],
            runtime_imported=payload["runtime_imported"],
            method_specific_soc_implemented=payload[
                "method_specific_soc_implemented"
            ],
            soc_derivatives_implemented=payload["soc_derivatives_implemented"],
            wavefunction_overlaps_implemented=payload[
                "wavefunction_overlaps_implemented"
            ],
            artifact_parser_validated=payload["artifact_parser_validated"],
            fresh_execution_observed=payload["fresh_execution_observed"],
        ).validate()


@dataclass(frozen=True)
class MolecularSOCAdmissionDossierV231:
    replay_dataset_fingerprint: str
    calculation_template_artifact: str
    environment_artifact: str
    artifacts: tuple[RawArtifactRecordV231, ...]
    receipts: tuple[CalculationReceiptV231, ...]
    trajectory_record_ids: tuple[str, ...]
    evidence: DerivedEvidenceBundleV231
    runtime_attestation: BackendRuntimeAttestationV231 | None = None

    def _maps(self, bundle_directory=None):
        artifacts = tuple(record.validate(bundle_directory) for record in self.artifacts)
        artifact_map = {record.name: record for record in artifacts}
        if len(artifact_map) != len(artifacts):
            raise ValueError("dossier artifact names must be unique.")
        if len({record.relative_path for record in artifacts}) != len(artifacts):
            raise ValueError("dossier artifact paths must be unique.")
        receipts = tuple(receipt.validate() for receipt in self.receipts)
        receipt_map = {receipt.record_id: receipt for receipt in receipts}
        if len(receipt_map) != len(receipts):
            raise ValueError("calculation receipt IDs must be unique.")
        return artifact_map, receipt_map

    def validate(self, *, bundle_directory=None, dataset=None, identity=None):
        _sha256_digest_v231(
            "replay_dataset_fingerprint", self.replay_dataset_fingerprint
        )
        self.evidence.validate()
        artifact_map, receipt_map = self._maps(bundle_directory)
        for name, role in (
            (self.calculation_template_artifact, "calculation_template"),
            (self.environment_artifact, "environment_lock"),
        ):
            if name not in artifact_map or artifact_map[name].role != role:
                raise ValueError(f"dossier lacks its {role} artifact.")
        input_names = set()
        output_names = set()
        for receipt in receipt_map.values():
            if receipt.input_artifact not in artifact_map:
                raise ValueError(f"receipt input artifact is missing: {receipt.record_id}.")
            if receipt.output_artifact not in artifact_map:
                raise ValueError(f"receipt output artifact is missing: {receipt.record_id}.")
            if artifact_map[receipt.input_artifact].role != "calculation_input":
                raise ValueError("receipt input must reference a calculation_input artifact.")
            if artifact_map[receipt.output_artifact].role != "calculation_output":
                raise ValueError("receipt output must reference a calculation_output artifact.")
            input_names.add(receipt.input_artifact)
            output_names.add(receipt.output_artifact)
        if len(input_names) != len(receipt_map) or len(output_names) != len(receipt_map):
            raise ValueError("each calculation receipt requires distinct input/output artifacts.")
        if len(self.trajectory_record_ids) != len(set(self.trajectory_record_ids)):
            raise ValueError("trajectory receipt IDs must be unique.")
        for record_id in self.trajectory_record_ids:
            if record_id not in receipt_map or receipt_map[record_id].role != "trajectory":
                raise ValueError("trajectory_record_ids must reference trajectory receipts.")

        evidence_artifacts = {
            self.evidence.reference.computed_artifact,
            *self.evidence.basis.source_artifacts,
            *self.evidence.method.source_artifacts,
            *self.evidence.frame.source_artifacts,
        }
        if not evidence_artifacts.issubset(output_names):
            raise ValueError("derived evidence must reference calculation output artifacts.")
        reference_name = self.evidence.reference.reference_artifact
        if (
            reference_name not in artifact_map
            or artifact_map[reference_name].role != "independent_reference"
        ):
            raise ValueError("independent reference artifact is absent or misclassified.")
        if reference_name in output_names:
            raise ValueError("independent reference cannot be a calculation output artifact.")

        output_to_receipt = {
            receipt.output_artifact: receipt for receipt in receipt_map.values()
        }
        for label, artifact in zip(
            self.evidence.basis.labels, self.evidence.basis.source_artifacts
        ):
            receipt = output_to_receipt[artifact]
            if receipt.role != "basis" or receipt.basis != label:
                raise ValueError("basis evidence labels disagree with calculation receipts.")
        for label, artifact in zip(
            self.evidence.method.labels, self.evidence.method.source_artifacts
        ):
            receipt = output_to_receipt[artifact]
            if receipt.role != "method" or receipt.electronic_method != label:
                raise ValueError("method evidence labels disagree with calculation receipts.")
        expected_frame_roles = ("frame_base", "frame_translation", "frame_rotation")
        actual_frame_roles = tuple(
            output_to_receipt[artifact].role
            for artifact in self.evidence.frame.source_artifacts
        )
        if actual_frame_roles != expected_frame_roles:
            raise ValueError("frame evidence artifacts have incorrect receipt roles.")

        if identity is not None:
            if not isinstance(identity, MolecularSOCBackendIdentityV230):
                raise TypeError("dossier identity must be MolecularSOCBackendIdentityV230.")
            identity.validate()
            template = artifact_map[self.calculation_template_artifact]
            environment = artifact_map[self.environment_artifact]
            if template.sha256 != identity.calculation_input_sha256:
                raise ValueError("calculation template hash disagrees with backend identity.")
            if environment.sha256 != identity.environment_sha256:
                raise ValueError("environment artifact hash disagrees with backend identity.")
            for receipt in receipt_map.values():
                common = (
                    receipt.backend_name == identity.backend_name
                    and receipt.backend_version == identity.backend_version
                    and receipt.source_kind == identity.source_kind
                    and receipt.soc_operator == identity.soc_operator
                    and receipt.derivative_method == identity.derivative_method
                )
                if not common:
                    raise ValueError("calculation receipt disagrees with backend identity.")
                if receipt.role not in {"basis", "method"} and (
                    receipt.electronic_method != identity.electronic_method
                    or receipt.basis != identity.basis
                ):
                    raise ValueError("non-ladder receipt changed method or basis.")

        if dataset is not None:
            if dataset.dataset_fingerprint != self.replay_dataset_fingerprint:
                raise ValueError("dossier is bound to a different replay dataset.")
            if len(self.trajectory_record_ids) != len(dataset.q):
                raise ValueError("trajectory receipts do not cover every replay record.")
            for index, record_id in enumerate(self.trajectory_record_ids):
                q = np.asarray(receipt_map[record_id].q_bohr, dtype=float)
                if q.shape != dataset.q[index].shape or not np.array_equal(q, dataset.q[index]):
                    raise ValueError("trajectory receipt coordinate disagrees with replay.")
            self.evidence.tracking.validate(
                nrecord=len(dataset.q), nstate=dataset.H_spin_free.shape[1]
            )

        if self.runtime_attestation is not None:
            attestation = self.runtime_attestation.validate()
            if (
                attestation.runtime_probe_artifact not in artifact_map
                or artifact_map[attestation.runtime_probe_artifact].role
                != "runtime_probe"
            ):
                raise ValueError("runtime attestation lacks its probe artifact.")
            if attestation.environment_sha256 != artifact_map[self.environment_artifact].sha256:
                raise ValueError("runtime attestation environment hash disagrees with dossier.")
        return self

    def derived_v230_evidence(self, dataset):
        self.validate(dataset=dataset)
        return self.evidence.derive_v230(dataset.overlaps)

    def as_dict(self):
        self.validate()
        return {
            "format": "gaussian-nadyn-molecular-soc-admission-dossier",
            "format_version": 1,
            "release": "v0.23.1",
            "replay_dataset_fingerprint": self.replay_dataset_fingerprint,
            "calculation_template_artifact": self.calculation_template_artifact,
            "environment_artifact": self.environment_artifact,
            "artifacts": [record.as_dict() for record in self.artifacts],
            "receipts": [receipt.as_dict() for receipt in self.receipts],
            "trajectory_record_ids": list(self.trajectory_record_ids),
            "evidence": self.evidence.as_dict(),
            "runtime_attestation": (
                None
                if self.runtime_attestation is None
                else self.runtime_attestation.as_dict()
            ),
        }

    def fingerprint(self):
        return hashlib.sha256(_canonical_json_bytes_v231(self.as_dict())).hexdigest()

    @classmethod
    def from_dict(cls, payload):
        if payload.get("format") != "gaussian-nadyn-molecular-soc-admission-dossier":
            raise ValueError("molecular SOC admission dossier format mismatch.")
        if payload.get("format_version") != 1:
            raise ValueError("molecular SOC admission dossier version mismatch.")
        if payload.get("release") != "v0.23.1":
            raise ValueError("molecular SOC admission dossier release mismatch.")
        runtime = payload.get("runtime_attestation")
        return cls(
            replay_dataset_fingerprint=payload["replay_dataset_fingerprint"],
            calculation_template_artifact=payload[
                "calculation_template_artifact"
            ],
            environment_artifact=payload["environment_artifact"],
            artifacts=tuple(
                RawArtifactRecordV231.from_dict(item) for item in payload["artifacts"]
            ),
            receipts=tuple(
                CalculationReceiptV231.from_dict(item) for item in payload["receipts"]
            ),
            trajectory_record_ids=tuple(payload["trajectory_record_ids"]),
            evidence=DerivedEvidenceBundleV231.from_dict(payload["evidence"]),
            runtime_attestation=(
                None
                if runtime is None
                else BackendRuntimeAttestationV231.from_dict(runtime)
            ),
        ).validate()


def write_molecular_soc_dossier_v231(
    bundle_directory,
    dossier,
    *,
    dataset=None,
    identity=None,
    overwrite=False,
):
    bundle_directory = Path(bundle_directory)
    bundle_directory.mkdir(parents=True, exist_ok=True)
    dossier.validate(
        bundle_directory=bundle_directory, dataset=dataset, identity=identity
    )
    path = bundle_directory / DOSSIER_NAME_V231
    if path.exists() and not overwrite:
        raise FileExistsError("molecular SOC admission dossier already exists.")
    payload = dossier.as_dict()
    manifest = {**payload, "dossier_fingerprint": dossier.fingerprint()}
    temporary = bundle_directory / f".{DOSSIER_NAME_V231}.tmp"
    temporary.write_bytes(_canonical_json_bytes_v231(manifest))
    os.replace(temporary, path)
    return load_molecular_soc_dossier_v231(
        path, dataset=dataset, identity=identity
    )


def load_molecular_soc_dossier_v231(path, *, dataset=None, identity=None):
    path = Path(path)
    if path.is_dir():
        path = path / DOSSIER_NAME_V231
    payload = json.loads(path.read_text(encoding="utf-8"))
    stored = payload.pop("dossier_fingerprint", None)
    dossier = MolecularSOCAdmissionDossierV231.from_dict(payload)
    if stored != dossier.fingerprint():
        raise ValueError("molecular SOC admission dossier integrity check failed.")
    return dossier.validate(
        bundle_directory=path.parent, dataset=dataset, identity=identity
    )
