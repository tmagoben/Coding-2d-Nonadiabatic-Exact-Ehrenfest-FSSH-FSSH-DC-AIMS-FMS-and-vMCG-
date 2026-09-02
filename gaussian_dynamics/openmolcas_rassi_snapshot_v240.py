"""Strict raw-artifact parser for the v0.24.0 OpenMolcas RASSI-SO intake.

Numerical arrays are read from a versioned exporter product, but that product is
cryptographically bound to the native input, text output, and ``rassi.h5`` file.
Protocol fixtures exercise this code path and remain distinguishable from external
executions at the type and admission layers.
"""

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path

import numpy as np

from .openmolcas_rassi_protocol_v240 import (
    OPENMOLCAS_ADAPTER_NAME_V240,
    OPENMOLCAS_ADAPTER_VERSION_V240,
    OPENMOLCAS_EXPORT_SCHEMA_V240,
    OPENMOLCAS_MANIFEST_SCHEMA_V240,
    OpenMolcasRASSIProtocolV240,
    openmolcas_protocol_from_dict_v240,
)


OPENMOLCAS_MANIFEST_NAME_V240 = "gnd_openmolcas_bundle_v240.json"
OPENMOLCAS_EXPORT_NAME_V240 = "gnd_rassi_export_v240.json"
OPENMOLCAS_INPUT_NAME_V240 = "openmolcas.input"
OPENMOLCAS_OUTPUT_NAME_V240 = "openmolcas.output"
OPENMOLCAS_HDF5_NAME_V240 = "rassi.h5"
OPENMOLCAS_VALIDATION_NAME_V240 = "gnd_external_validation_v240.json"
OPENMOLCAS_VALIDATION_ARTIFACT_DIRECTORY_V240 = "validation_artifacts"
PROTOCOL_FIXTURE_MARKER_V240 = "GND-V0240-PROTOCOL-FIXTURE-NO-OPENMOLCAS-EXECUTION"
HDF5_MAGIC_V240 = b"\x89HDF\r\n\x1a\n"
# The bundle parser binds exported arrays to native files by digest, but v0.24.0
# does not yet independently reconstruct the SOC matrices from OpenMolcas-native
# HDF5/text datasets.  Admission must remain closed until that cross-parser exists.
NATIVE_OPENMOLCAS_NUMERIC_CROSSCHECK_V240 = False

_CONVERGENCE_KEYS_V240 = {
    "gateway",
    "seward",
    "scf",
    "rasscf_singlet",
    "rasscf_triplet",
    "caspt2_singlet",
    "caspt2_triplet",
    "rassi_so",
}


def sha256_file_v240(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _sha256_v240(name, value):
    text = str(value)
    if len(text) != 64 or any(item not in "0123456789abcdef" for item in text):
        raise ValueError(f"{name} must be a lowercase SHA-256 digest.")
    return text


def _strict_json_v240(path):
    return json.loads(
        Path(path).read_text(encoding="utf-8"),
        parse_constant=lambda value: (_ for _ in ()).throw(
            ValueError(f"non-finite JSON constant {value!r} is forbidden.")
        ),
    )


def _complex_matrix_v240(payload, real_name, imag_name, shape):
    real = np.asarray(payload[real_name], dtype=float)
    imag = np.asarray(payload[imag_name], dtype=float)
    if real.shape != shape or imag.shape != shape:
        raise ValueError(f"{real_name}/{imag_name} must have shape {shape}.")
    value = real + 1j * imag
    if not np.all(np.isfinite(value)):
        raise ValueError(f"{real_name}/{imag_name} contain non-finite data.")
    return value


@dataclass(frozen=True)
class OpenMolcasArtifactRecordV240:
    record_id: str
    relative_directory: str
    input_sha256: str
    output_sha256: str
    rassi_h5_sha256: str
    export_sha256: str

    def validate(self):
        if not self.record_id or not self.relative_directory:
            raise ValueError("artifact record identifiers cannot be empty.")
        path = Path(self.relative_directory)
        if path.is_absolute() or ".." in path.parts or len(path.parts) != 2:
            raise ValueError("record directory must be records/<record_id>.")
        if path.parts != ("records", self.record_id):
            raise ValueError("record directory and record_id disagree.")
        for name in (
            "input_sha256",
            "output_sha256",
            "rassi_h5_sha256",
            "export_sha256",
        ):
            _sha256_v240(name, getattr(self, name))
        return self

    @classmethod
    def from_dict(cls, payload):
        if not isinstance(payload, dict) or set(payload) != set(cls.__dataclass_fields__):
            raise ValueError("OpenMolcas artifact record fields are incomplete or unknown.")
        return cls(**payload).validate()


@dataclass(frozen=True)
class OpenMolcasBundleManifestV240:
    schema: str
    source_kind: str
    protocol: OpenMolcasRASSIProtocolV240
    adapter_name: str
    adapter_version: str
    exporter_name: str
    exporter_version: str
    environment_sha256: str
    validation_artifact: str
    validation_sha256: str
    records: tuple[OpenMolcasArtifactRecordV240, ...]

    def validate(self):
        if self.schema != OPENMOLCAS_MANIFEST_SCHEMA_V240:
            raise ValueError("OpenMolcas bundle manifest schema mismatch.")
        if self.source_kind not in {
            "protocol_fixture",
            "external_ab_initio_snapshot",
        }:
            raise ValueError("unsupported OpenMolcas bundle source_kind.")
        self.protocol.validate()
        if self.adapter_name != OPENMOLCAS_ADAPTER_NAME_V240:
            raise ValueError("OpenMolcas adapter name mismatch.")
        if self.adapter_version != OPENMOLCAS_ADAPTER_VERSION_V240:
            raise ValueError("OpenMolcas adapter version mismatch.")
        if not self.exporter_name.strip() or not self.exporter_version.strip():
            raise ValueError("exporter identity cannot be empty.")
        _sha256_v240("environment_sha256", self.environment_sha256)
        if self.validation_artifact != OPENMOLCAS_VALIDATION_NAME_V240:
            raise ValueError("external validation artifact name mismatch.")
        _sha256_v240("validation_sha256", self.validation_sha256)
        records = tuple(item.validate() for item in self.records)
        identifiers = tuple(item.record_id for item in records)
        if identifiers != self.protocol.expected_record_ids():
            raise ValueError("bundle record inventory/order is not the frozen protocol.")
        return self

    @classmethod
    def from_dict(cls, payload):
        expected = {
            "schema",
            "source_kind",
            "protocol",
            "adapter_name",
            "adapter_version",
            "exporter_name",
            "exporter_version",
            "environment_sha256",
            "validation_artifact",
            "validation_sha256",
            "records",
        }
        if not isinstance(payload, dict) or set(payload) != expected:
            raise ValueError("OpenMolcas bundle manifest fields are incomplete or unknown.")
        return cls(
            schema=payload["schema"],
            source_kind=payload["source_kind"],
            protocol=openmolcas_protocol_from_dict_v240(payload["protocol"]),
            adapter_name=payload["adapter_name"],
            adapter_version=payload["adapter_version"],
            exporter_name=payload["exporter_name"],
            exporter_version=payload["exporter_version"],
            environment_sha256=payload["environment_sha256"],
            validation_artifact=payload["validation_artifact"],
            validation_sha256=payload["validation_sha256"],
            records=tuple(
                OpenMolcasArtifactRecordV240.from_dict(item)
                for item in payload["records"]
            ),
        ).validate()


@dataclass(frozen=True)
class ParsedOpenMolcasRecordV240:
    record_id: str
    geometry_bohr: np.ndarray
    H_spin_free: np.ndarray
    H_soc: np.ndarray
    reference_overlap: np.ndarray
    convergence: dict
    artifact_digests: dict
    native_hdf5: bool
    native_output_completed: bool

    @property
    def H_total(self):
        return self.H_spin_free + self.H_soc


@dataclass(frozen=True)
class ParsedOpenMolcasBundleV240:
    directory: Path
    manifest: OpenMolcasBundleManifestV240
    manifest_sha256: str
    records: tuple[ParsedOpenMolcasRecordV240, ...]
    exact_artifact_inventory: bool
    parser_executed: bool

    @property
    def source_kind(self):
        return self.manifest.source_kind

    @property
    def protocol(self):
        return self.manifest.protocol

    @property
    def record_map(self):
        return {item.record_id: item for item in self.records}

    @property
    def native_openmolcas_execution(self):
        return bool(
            self.source_kind == "external_ab_initio_snapshot"
            and all(item.native_hdf5 for item in self.records)
            and all(item.native_output_completed for item in self.records)
        )

    @property
    def fingerprint(self):
        digest = hashlib.sha256()
        digest.update(self.manifest_sha256.encode("ascii"))
        digest.update(self.manifest.validation_sha256.encode("ascii"))
        for record in self.records:
            for name in sorted(record.artifact_digests):
                digest.update(name.encode("utf-8"))
                digest.update(record.artifact_digests[name].encode("ascii"))
        return digest.hexdigest()


class OpenMolcasRASSISnapshotParserV240:
    """Exact-type parser used by the v0.24.0 admission trust policy."""

    parser_name = OPENMOLCAS_ADAPTER_NAME_V240
    parser_version = OPENMOLCAS_ADAPTER_VERSION_V240

    @staticmethod
    def _validate_input(text, protocol):
        upper = text.upper()
        required = (
            "&GATEWAY",
            "&SEWARD",
            "AMFI",
            "&RASSCF",
            "&CASPT2",
            "&RASSI",
            "NROFJOBIPH",
            "SPINORBIT",
            "EJOB",
            protocol.basis.upper(),
            "UNIT=BOHR",
            "GROUP=NOSYM",
            "R02O",
            "INACTIVE=1",
            "CIROOT=1 1 1",
            ">> COPY",
        )
        if any(marker not in upper for marker in required):
            raise ValueError("OpenMolcas input lacks a frozen protocol module/keyword.")
        lines = [line.strip() for line in text.splitlines()]
        begin_marker = "* GND_GEOMETRY_BOHR_BEGIN"
        end_marker = "* GND_GEOMETRY_BOHR_END"
        if lines.count(begin_marker) != 1 or lines.count(end_marker) != 1:
            raise ValueError("OpenMolcas input lacks unique geometry binding markers.")
        begin = lines.index(begin_marker)
        end = lines.index(end_marker)
        if end <= begin + 1:
            raise ValueError("OpenMolcas input geometry binding is empty.")
        entries = lines[begin + 1 : end]
        if len(entries) != len(protocol.atom_symbols):
            raise ValueError("OpenMolcas input geometry has the wrong atom count.")
        symbols = []
        coordinates = []
        for entry in entries:
            fields = entry.split()
            if len(fields) != 4:
                raise ValueError("OpenMolcas input geometry row is malformed.")
            symbols.append(fields[0])
            coordinates.append([float(value) for value in fields[1:]])
        if tuple(symbols) != protocol.atom_symbols:
            raise ValueError("OpenMolcas input nuclear order differs from the protocol.")
        geometry = np.asarray(coordinates, dtype=float)
        if not np.all(np.isfinite(geometry)):
            raise ValueError("OpenMolcas input geometry contains non-finite data.")
        return geometry

    @staticmethod
    def _output_status(text, source_kind, protocol):
        if source_kind == "protocol_fixture":
            if PROTOCOL_FIXTURE_MARKER_V240 not in text:
                raise ValueError("protocol fixture lacks its mandatory non-execution marker.")
            return False
        if PROTOCOL_FIXTURE_MARKER_V240 in text:
            raise ValueError("protocol fixture cannot be relabeled as external evidence.")
        upper = text.upper()
        required = (
            "OPENMOLCAS",
            protocol.backend_version.upper(),
            "RASSI",
            "SPIN-ORBIT",
            "HAPPY LANDING",
        )
        if any(marker not in upper for marker in required):
            raise ValueError("external OpenMolcas output lacks required completion markers.")
        return True

    def parse_bundle_v240(self, directory):
        root = Path(directory).resolve()
        if not root.is_dir():
            raise ValueError("OpenMolcas bundle directory does not exist.")
        manifest_path = root / OPENMOLCAS_MANIFEST_NAME_V240
        if manifest_path.is_symlink() or not manifest_path.is_file():
            raise ValueError("OpenMolcas bundle requires a regular manifest file.")
        manifest = OpenMolcasBundleManifestV240.from_dict(
            _strict_json_v240(manifest_path)
        )
        expected_files = {manifest_path.resolve()}
        validation_path = root / manifest.validation_artifact
        if validation_path.is_symlink() or not validation_path.is_file():
            raise ValueError("bundle lacks its regular external-validation artifact.")
        if sha256_file_v240(validation_path) != manifest.validation_sha256:
            raise ValueError("external-validation artifact digest mismatch.")
        expected_files.add(validation_path.resolve())
        validation_payload = _strict_json_v240(validation_path)
        validation_digest_fields = (
            "basis_artifact_sha256",
            "method_artifact_sha256",
            "frame_artifact_sha256",
            "tracking_artifact_sha256",
        )
        if "independent_artifact_sha256" not in validation_payload or any(
            name not in validation_payload for name in validation_digest_fields
        ):
            raise ValueError("external validation lacks its raw-artifact inventory.")
        validation_digests = [validation_payload["independent_artifact_sha256"]]
        for name in validation_digest_fields:
            values = validation_payload[name]
            if not isinstance(values, list):
                raise ValueError("external validation raw-artifact inventory is malformed.")
            validation_digests.extend(values)
        if len(validation_digests) != len(set(validation_digests)):
            raise ValueError("external validation raw artifacts must be unique.")
        validation_artifact_directory = (
            root / OPENMOLCAS_VALIDATION_ARTIFACT_DIRECTORY_V240
        )
        if (
            validation_artifact_directory.is_symlink()
            or not validation_artifact_directory.is_dir()
        ):
            raise ValueError("external validation raw-artifact directory is absent.")
        for digest in validation_digests:
            _sha256_v240("external validation artifact sha256", digest)
            artifact_path = validation_artifact_directory / f"{digest}.artifact"
            if artifact_path.is_symlink() or not artifact_path.is_file():
                raise ValueError("external validation raw artifact is absent or a symlink.")
            if sha256_file_v240(artifact_path) != digest:
                raise ValueError("external validation raw artifact digest mismatch.")
            expected_files.add(artifact_path.resolve())
        parsed = []
        for artifact in manifest.records:
            record_directory = root / artifact.relative_directory
            if record_directory.is_symlink() or not record_directory.is_dir():
                raise ValueError("record directory must be a regular in-bundle directory.")
            paths = {
                "input": record_directory / OPENMOLCAS_INPUT_NAME_V240,
                "output": record_directory / OPENMOLCAS_OUTPUT_NAME_V240,
                "rassi_h5": record_directory / OPENMOLCAS_HDF5_NAME_V240,
                "export": record_directory / OPENMOLCAS_EXPORT_NAME_V240,
            }
            if any(path.is_symlink() or not path.is_file() for path in paths.values()):
                raise ValueError("record artifact is absent, non-regular, or a symlink.")
            expected_files.update(path.resolve() for path in paths.values())
            digests = {name: sha256_file_v240(path) for name, path in paths.items()}
            declared = {
                "input": artifact.input_sha256,
                "output": artifact.output_sha256,
                "rassi_h5": artifact.rassi_h5_sha256,
                "export": artifact.export_sha256,
            }
            if digests != declared:
                raise ValueError("record artifact digest mismatch.")
            input_text = paths["input"].read_text(encoding="utf-8")
            output_text = paths["output"].read_text(encoding="utf-8")
            input_geometry = self._validate_input(input_text, manifest.protocol)
            completed = self._output_status(
                output_text, manifest.source_kind, manifest.protocol
            )
            native_hdf5 = paths["rassi_h5"].read_bytes()[:8] == HDF5_MAGIC_V240
            if manifest.source_kind == "external_ab_initio_snapshot" and not native_hdf5:
                raise ValueError("external snapshot rassi.h5 lacks the HDF5 file signature.")
            export = _strict_json_v240(paths["export"])
            expected_export = {
                "schema",
                "record_id",
                "protocol_fingerprint",
                "input_sha256",
                "output_sha256",
                "rassi_h5_sha256",
                "geometry_bohr",
                "state_labels",
                "H_spin_free_real",
                "H_spin_free_imag",
                "H_soc_real",
                "H_soc_imag",
                "reference_overlap_real",
                "reference_overlap_imag",
                "convergence",
            }
            if not isinstance(export, dict) or set(export) != expected_export:
                raise ValueError("OpenMolcas export fields are incomplete or unknown.")
            if export["schema"] != OPENMOLCAS_EXPORT_SCHEMA_V240:
                raise ValueError("OpenMolcas export schema mismatch.")
            if export["record_id"] != artifact.record_id:
                raise ValueError("export record_id disagrees with its directory.")
            if export["protocol_fingerprint"] != manifest.protocol.fingerprint():
                raise ValueError("export protocol fingerprint mismatch.")
            for key, export_name in (
                ("input", "input_sha256"),
                ("output", "output_sha256"),
                ("rassi_h5", "rassi_h5_sha256"),
            ):
                if export[export_name] != digests[key]:
                    raise ValueError("export is not bound to its native artifact set.")
            if tuple(export["state_labels"]) != manifest.protocol.state_order:
                raise ValueError("exported electronic-state order is not the protocol order.")
            geometry = np.asarray(export["geometry_bohr"], dtype=float)
            natom = len(manifest.protocol.atom_symbols)
            if geometry.shape != (natom, 3) or not np.all(np.isfinite(geometry)):
                raise ValueError("exported geometry has invalid shape or values.")
            if np.max(np.abs(geometry - input_geometry)) > 1.0e-12:
                raise ValueError("exported geometry is not bound to the native input.")
            nstate = len(manifest.protocol.state_order)
            H0 = _complex_matrix_v240(
                export, "H_spin_free_real", "H_spin_free_imag", (nstate, nstate)
            )
            Hso = _complex_matrix_v240(
                export, "H_soc_real", "H_soc_imag", (nstate, nstate)
            )
            overlap = _complex_matrix_v240(
                export,
                "reference_overlap_real",
                "reference_overlap_imag",
                (nstate, nstate),
            )
            convergence = export["convergence"]
            if not isinstance(convergence, dict) or set(convergence) != _CONVERGENCE_KEYS_V240:
                raise ValueError("exported convergence inventory is incomplete or unknown.")
            if any(type(value) is not bool for value in convergence.values()):
                raise TypeError("all convergence values must be native Booleans.")
            parsed.append(
                ParsedOpenMolcasRecordV240(
                    record_id=artifact.record_id,
                    geometry_bohr=geometry,
                    H_spin_free=H0,
                    H_soc=Hso,
                    reference_overlap=overlap,
                    convergence=dict(convergence),
                    artifact_digests=digests,
                    native_hdf5=bool(native_hdf5),
                    native_output_completed=bool(completed),
                )
            )
        observed_files = {
            path.resolve()
            for path in root.rglob("*")
            if path.is_file() or path.is_symlink()
        }
        if observed_files != expected_files:
            raise ValueError("OpenMolcas bundle contains an unknown or missing artifact.")
        return ParsedOpenMolcasBundleV240(
            directory=root,
            manifest=manifest,
            manifest_sha256=sha256_file_v240(manifest_path),
            records=tuple(parsed),
            exact_artifact_inventory=True,
            parser_executed=True,
        )
