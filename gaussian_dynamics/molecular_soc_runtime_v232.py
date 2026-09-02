"""Fail-closed runtime evidence types for molecular SOC admission in v0.23.2.

This module does not implement a spin--orbit method.  It defines the evidence that a
method-specific implementation must emit and the out-of-band trust policy against
which that evidence is checked.  In particular, an engine cannot make itself trusted
by declaring its own adapter name or capability flags.
"""

from dataclasses import asdict, dataclass
import json
from pathlib import Path
from typing import Mapping

import numpy as np

from .molecular_soc_contract_v230 import MolecularSOCBackendIdentityV230


CONVERGENCE_VOCABULARY_V232 = "gaussian-nadyn/molecular-soc-convergence/v2"
CONVERGENCE_METADATA_KEY_V232 = "v232_molecular_soc_convergence"
CONVERGENCE_STAGES_V232 = (
    "scf",
    "correlated_wavefunction",
    "state_interaction_soc",
    "spin_free_gradients",
    "soc_derivatives",
    "derivative_connections",
    "many_electron_overlaps",
)
RUNTIME_PROBE_FORMAT_V232 = "gaussian-nadyn-molecular-soc-runtime-probe"


def _required_text_v232(name, value):
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{name} must be a nonempty string.")
    return value


def _native_bool_v232(name, value):
    if not isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{name} must be Boolean.")
    return bool(value)


def _sha256_v232(name, value):
    if not isinstance(value, str) or len(value) != 64 or any(
        character not in "0123456789abcdef" for character in value
    ):
        raise ValueError(f"{name} must be a lowercase SHA-256 digest.")
    return value


@dataclass(frozen=True)
class ConvergenceMetadataV232:
    """Canonical, component-resolved convergence evidence.

    The v0.23.1 ``derivatives_converged`` spelling was intentionally not accepted:
    it cannot distinguish spin-free gradients, physical SOC derivatives, and
    derivative connections.  Adapters must translate native backend vocabulary into
    these seven explicit stages.
    """

    scf: bool
    correlated_wavefunction: bool
    state_interaction_soc: bool
    spin_free_gradients: bool
    soc_derivatives: bool
    derivative_connections: bool
    many_electron_overlaps: bool
    vocabulary: str = CONVERGENCE_VOCABULARY_V232

    def validate(self):
        if self.vocabulary != CONVERGENCE_VOCABULARY_V232:
            raise ValueError("unsupported molecular SOC convergence vocabulary.")
        for name in CONVERGENCE_STAGES_V232:
            _native_bool_v232(name, getattr(self, name))
        return self

    @property
    def complete(self):
        self.validate()
        return bool(all(getattr(self, name) for name in CONVERGENCE_STAGES_V232))

    def as_dict(self):
        self.validate()
        return {
            "vocabulary": self.vocabulary,
            "stages": {
                name: bool(getattr(self, name)) for name in CONVERGENCE_STAGES_V232
            },
        }

    @classmethod
    def from_mapping(cls, payload):
        if not isinstance(payload, Mapping):
            raise TypeError("convergence metadata must be a mapping.")
        if set(payload) != {"vocabulary", "stages"}:
            raise ValueError("convergence metadata has missing or unknown fields.")
        stages = payload["stages"]
        if not isinstance(stages, Mapping) or set(stages) != set(
            CONVERGENCE_STAGES_V232
        ):
            raise ValueError("convergence stage vocabulary is incomplete or ambiguous.")
        return cls(
            vocabulary=payload["vocabulary"],
            **{name: stages[name] for name in CONVERGENCE_STAGES_V232},
        ).validate()


def convergence_from_snapshot_v232(snapshot):
    """Read one unambiguous convergence namespace from a validated snapshot."""
    point_metadata = getattr(getattr(snapshot, "point", None), "metadata", None)
    snapshot_metadata = getattr(snapshot, "metadata", None)
    if not isinstance(point_metadata, Mapping) or not isinstance(
        snapshot_metadata, Mapping
    ):
        raise TypeError("snapshot and point metadata must be mappings.")
    if CONVERGENCE_METADATA_KEY_V232 in point_metadata:
        raise ValueError(
            "convergence namespace must occur only in snapshot metadata; duplicate "
            "point/snapshot declarations are forbidden."
        )
    if CONVERGENCE_METADATA_KEY_V232 not in snapshot_metadata:
        raise ValueError("snapshot lacks canonical v0.23.2 convergence metadata.")
    return ConvergenceMetadataV232.from_mapping(
        snapshot_metadata[CONVERGENCE_METADATA_KEY_V232]
    )


@dataclass(frozen=True)
class BackendMethodIdentityV232:
    """Exact backend, adapter, and electronic-method identity."""

    backend_name: str
    backend_version: str
    source_kind: str
    adapter_name: str
    adapter_version: str
    electronic_method: str
    basis: str
    active_space: str
    soc_operator: str
    scalar_relativistic_method: str
    derivative_method: str
    nac_convention: str

    def validate(self):
        for name in (
            "backend_name",
            "backend_version",
            "source_kind",
            "adapter_name",
            "adapter_version",
            "electronic_method",
            "basis",
            "active_space",
            "soc_operator",
            "scalar_relativistic_method",
            "derivative_method",
            "nac_convention",
        ):
            _required_text_v232(name, getattr(self, name))
        if self.source_kind not in {
            "external_ab_initio_snapshot",
            "live_ab_initio",
        }:
            raise ValueError("runtime identity requires an external or live source.")
        return self

    def matches_backend_identity(self, identity):
        self.validate()
        if not isinstance(identity, MolecularSOCBackendIdentityV230):
            raise TypeError("backend identity must be MolecularSOCBackendIdentityV230.")
        identity.validate()
        return bool(
            self.backend_name == identity.backend_name
            and self.backend_version == identity.backend_version
            and self.source_kind == identity.source_kind
            and self.electronic_method == identity.electronic_method
            and self.basis == identity.basis
            and self.active_space == identity.active_space
            and self.soc_operator == identity.soc_operator
            and self.scalar_relativistic_method
            == identity.scalar_relativistic_method
            and self.derivative_method == identity.derivative_method
        )

    def as_dict(self):
        self.validate()
        return asdict(self)

    @classmethod
    def from_dict(cls, payload):
        if not isinstance(payload, Mapping) or set(payload) != {
            "backend_name",
            "backend_version",
            "source_kind",
            "adapter_name",
            "adapter_version",
            "electronic_method",
            "basis",
            "active_space",
            "soc_operator",
            "scalar_relativistic_method",
            "derivative_method",
            "nac_convention",
        }:
            raise ValueError("backend method identity fields are incomplete or unknown.")
        return cls(**dict(payload)).validate()


@dataclass(frozen=True)
class RuntimeProbeRecordV232:
    method_identity: BackendMethodIdentityV232
    environment_sha256: str
    calculation_input_sha256: str
    replay_dataset_fingerprint: str
    runtime_imported: bool
    format: str = RUNTIME_PROBE_FORMAT_V232
    format_version: int = 2

    def validate(self):
        if self.format != RUNTIME_PROBE_FORMAT_V232 or self.format_version != 2:
            raise ValueError("runtime probe format mismatch.")
        self.method_identity.validate()
        _sha256_v232("environment_sha256", self.environment_sha256)
        _sha256_v232("calculation_input_sha256", self.calculation_input_sha256)
        _sha256_v232(
            "replay_dataset_fingerprint", self.replay_dataset_fingerprint
        )
        _native_bool_v232("runtime_imported", self.runtime_imported)
        return self

    def as_dict(self):
        self.validate()
        return {
            "format": self.format,
            "format_version": self.format_version,
            "method_identity": self.method_identity.as_dict(),
            "environment_sha256": self.environment_sha256,
            "calculation_input_sha256": self.calculation_input_sha256,
            "replay_dataset_fingerprint": self.replay_dataset_fingerprint,
            "runtime_imported": bool(self.runtime_imported),
        }

    @classmethod
    def from_dict(cls, payload):
        if not isinstance(payload, Mapping) or set(payload) != {
            "format",
            "format_version",
            "method_identity",
            "environment_sha256",
            "calculation_input_sha256",
            "replay_dataset_fingerprint",
            "runtime_imported",
        }:
            raise ValueError("runtime probe fields are incomplete or unknown.")
        return cls(
            format=payload["format"],
            format_version=payload["format_version"],
            method_identity=BackendMethodIdentityV232.from_dict(
                payload["method_identity"]
            ),
            environment_sha256=payload["environment_sha256"],
            calculation_input_sha256=payload["calculation_input_sha256"],
            replay_dataset_fingerprint=payload["replay_dataset_fingerprint"],
            runtime_imported=payload["runtime_imported"],
        ).validate()


def load_runtime_probe_v232(path):
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    return RuntimeProbeRecordV232.from_dict(payload)


@dataclass(frozen=True)
class ReceiptExecutionEvidenceV232:
    record_id: str
    output_artifact: str
    convergence: ConvergenceMetadataV232

    def validate(self):
        _required_text_v232("record_id", self.record_id)
        _required_text_v232("output_artifact", self.output_artifact)
        self.convergence.validate()
        return self


@dataclass(frozen=True)
class BackendArtifactValidationProofV232:
    """Typed result of executing a trusted raw-artifact parser."""

    method_identity: BackendMethodIdentityV232
    parser_name: str
    parser_version: str
    execution_challenge: str
    replay_dataset_fingerprint: str
    dossier_fingerprint: str
    environment_sha256: str
    calculation_input_sha256: str
    runtime_probe_artifact: str
    runtime_probe_sha256: str
    parsed_output_artifacts: tuple[tuple[str, str], ...]
    receipt_evidence: tuple[ReceiptExecutionEvidenceV232, ...]
    parser_executed: bool
    fresh_execution_observed: bool

    def validate(self):
        self.method_identity.validate()
        _required_text_v232("parser_name", self.parser_name)
        _required_text_v232("parser_version", self.parser_version)
        for name in (
            "execution_challenge",
            "replay_dataset_fingerprint",
            "dossier_fingerprint",
            "environment_sha256",
            "calculation_input_sha256",
            "runtime_probe_sha256",
        ):
            _sha256_v232(name, getattr(self, name))
        _required_text_v232("runtime_probe_artifact", self.runtime_probe_artifact)
        artifact_names = []
        for name, digest in self.parsed_output_artifacts:
            artifact_names.append(_required_text_v232("parsed artifact name", name))
            _sha256_v232("parsed artifact sha256", digest)
        if len(artifact_names) != len(set(artifact_names)):
            raise ValueError("parsed output artifacts must be unique.")
        receipts = tuple(item.validate() for item in self.receipt_evidence)
        if len({item.record_id for item in receipts}) != len(receipts):
            raise ValueError("receipt execution evidence IDs must be unique.")
        _native_bool_v232("parser_executed", self.parser_executed)
        _native_bool_v232(
            "fresh_execution_observed", self.fresh_execution_observed
        )
        return self


@dataclass(frozen=True)
class BackendAdmissionPolicyV232:
    """Out-of-band trust anchor; it must not be read from the dossier or engine."""

    expected_identity: BackendMethodIdentityV232
    trusted_validator_type: type
    parser_name: str
    parser_version: str

    def validate(self):
        self.expected_identity.validate()
        if not isinstance(self.trusted_validator_type, type):
            raise TypeError("trusted_validator_type must be a concrete Python type.")
        _required_text_v232("parser_name", self.parser_name)
        _required_text_v232("parser_version", self.parser_version)
        return self
