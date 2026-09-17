"""v0.29.4 pre-ML data and evaluation contract for Gaussian spawning.

This module deliberately does *not* contain an ML model.  Its job is to freeze the
interface between a deterministic physics spawning oracle and any future surrogate.
The surrogate may rank/propose candidates; the physics oracle remains authoritative.
"""

from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
import json
import math
from typing import Any, Iterable, Mapping, Sequence

import numpy as np

SCHEMA_VERSION = "0.29.4"
_ALLOWED_SPLITS = frozenset({"train", "validation", "test", "audit"})
_ALLOWED_ACTIONS = frozenset({"spawn", "no_spawn"})

# Future models are evaluated on physics-relevant failure modes, not accuracy alone.
REQUIRED_EVALUATION_METRICS = (
    "precision",
    "recall",
    "false_negative_rate",
    "brier_score",
    "mean_score_error",
    "event_time_mae",
)


def _finite(value: float, name: str) -> float:
    value = float(value)
    if not math.isfinite(value):
        raise ValueError(f"{name} must be finite")
    return value


def _finite_tuple(values: Sequence[float], name: str) -> tuple[float, ...]:
    return tuple(_finite(value, f"{name}[{i}]") for i, value in enumerate(values))


def _nonempty(text: str, name: str) -> str:
    if not isinstance(text, str) or not text.strip():
        raise ValueError(f"{name} must be a non-empty string")
    return text.strip()


def canonical_json(payload: Mapping[str, Any] | Sequence[Any]) -> str:
    """Return a byte-stable JSON representation used by all v0.29.4 hashes."""
    return json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    )


def fingerprint(payload: Mapping[str, Any] | Sequence[Any]) -> str:
    """SHA-256 of :func:`canonical_json`."""
    return sha256(canonical_json(payload).encode("utf-8")).hexdigest()


@dataclass(frozen=True)
class EventIdentity:
    """Stable identity and leakage group for one candidate evaluation."""

    trajectory_id: str
    family_id: str
    event_id: str
    packet_id: str
    step: int
    time_au: float
    parent_packet_id: str | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "trajectory_id", _nonempty(self.trajectory_id, "trajectory_id"))
        object.__setattr__(self, "family_id", _nonempty(self.family_id, "family_id"))
        object.__setattr__(self, "event_id", _nonempty(self.event_id, "event_id"))
        object.__setattr__(self, "packet_id", _nonempty(self.packet_id, "packet_id"))
        if self.parent_packet_id is not None:
            object.__setattr__(
                self,
                "parent_packet_id",
                _nonempty(self.parent_packet_id, "parent_packet_id"),
            )
        if int(self.step) != self.step or self.step < 0:
            raise ValueError("step must be a non-negative integer")
        object.__setattr__(self, "step", int(self.step))
        time_au = _finite(self.time_au, "time_au")
        if time_au < 0.0:
            raise ValueError("time_au must be non-negative")
        object.__setattr__(self, "time_au", time_au)

    @property
    def leakage_group(self) -> tuple[str, str]:
        """Family + trajectory is the minimum allowed split granularity."""
        return (self.family_id, self.trajectory_id)

    def to_payload(self) -> dict[str, Any]:
        return {
            "trajectory_id": self.trajectory_id,
            "family_id": self.family_id,
            "event_id": self.event_id,
            "packet_id": self.packet_id,
            "parent_packet_id": self.parent_packet_id,
            "step": self.step,
            "time_au": self.time_au,
        }


@dataclass(frozen=True)
class InvariantFeatures:
    """Gauge/packet-frame invariant scalars presented to a future surrogate.

    No raw adiabatic eigenvectors, NAC vectors, SOC matrix entries, phases, or gauge
    labels are admitted here.  Quantities are intentionally reduced to invariant
    populations, spectra, norms, singular values, and conditioning diagnostics.
    """

    populations: tuple[float, ...]
    energy_gaps: tuple[float, ...] = ()
    force_norms: tuple[float, ...] = ()
    nac_frobenius: float | None = None
    soc_singular_values: tuple[float, ...] = ()
    berry_curvature_frobenius: float | None = None
    gram_min_eigenvalue: float | None = None
    gram_condition_number: float | None = None
    max_packet_overlap_abs: float | None = None
    nuclear_speed: float | None = None
    kinetic_energy: float | None = None
    residual_norm: float | None = None
    novelty: float | None = None
    projection_loss_estimate: float | None = None
    dt_au: float | None = None

    def __post_init__(self) -> None:
        populations = _finite_tuple(self.populations, "populations")
        if len(populations) < 2:
            raise ValueError("populations must contain at least two electronic states")
        if any(value < -1e-12 for value in populations):
            raise ValueError("populations must be non-negative")
        if not math.isclose(sum(populations), 1.0, rel_tol=0.0, abs_tol=1e-9):
            raise ValueError("populations must sum to 1 within 1e-9")
        object.__setattr__(self, "populations", populations)

        object.__setattr__(self, "energy_gaps", _finite_tuple(self.energy_gaps, "energy_gaps"))
        object.__setattr__(self, "force_norms", _finite_tuple(self.force_norms, "force_norms"))
        object.__setattr__(
            self,
            "soc_singular_values",
            _finite_tuple(self.soc_singular_values, "soc_singular_values"),
        )
        if any(value < -1e-14 for value in self.energy_gaps):
            raise ValueError("energy_gaps must be non-negative ordered gaps")
        if any(value < -1e-14 for value in self.force_norms):
            raise ValueError("force_norms must be non-negative")
        if any(value < -1e-14 for value in self.soc_singular_values):
            raise ValueError("soc_singular_values must be non-negative")

        nonnegative = {
            "nac_frobenius": self.nac_frobenius,
            "berry_curvature_frobenius": self.berry_curvature_frobenius,
            "gram_min_eigenvalue": self.gram_min_eigenvalue,
            "gram_condition_number": self.gram_condition_number,
            "max_packet_overlap_abs": self.max_packet_overlap_abs,
            "nuclear_speed": self.nuclear_speed,
            "kinetic_energy": self.kinetic_energy,
            "residual_norm": self.residual_norm,
            "projection_loss_estimate": self.projection_loss_estimate,
            "dt_au": self.dt_au,
        }
        for name, value in nonnegative.items():
            if value is None:
                continue
            checked = _finite(value, name)
            if checked < -1e-14:
                raise ValueError(f"{name} must be non-negative")
            object.__setattr__(self, name, checked)

        if self.novelty is not None:
            novelty = _finite(self.novelty, "novelty")
            if not 0.0 <= novelty <= 1.0:
                raise ValueError("novelty must lie in [0, 1]")
            object.__setattr__(self, "novelty", novelty)
        if self.max_packet_overlap_abs is not None and self.max_packet_overlap_abs > 1.0 + 1e-12:
            raise ValueError("max_packet_overlap_abs cannot exceed 1")
        if self.gram_condition_number is not None and self.gram_condition_number < 1.0 - 1e-12:
            raise ValueError("gram_condition_number must be >= 1")

    def to_payload(self) -> dict[str, Any]:
        payload: dict[str, Any] = {
            "populations": list(self.populations),
            "energy_gaps": list(self.energy_gaps),
            "force_norms": list(self.force_norms),
            "soc_singular_values": list(self.soc_singular_values),
        }
        for name in (
            "nac_frobenius",
            "berry_curvature_frobenius",
            "gram_min_eigenvalue",
            "gram_condition_number",
            "max_packet_overlap_abs",
            "nuclear_speed",
            "kinetic_energy",
            "residual_norm",
            "novelty",
            "projection_loss_estimate",
            "dt_au",
        ):
            value = getattr(self, name)
            if value is not None:
                payload[name] = value
        return payload


@dataclass(frozen=True)
class OracleDecision:
    """Exact physics-oracle outcome used as supervision."""

    action: str
    residual_gain: float
    normalized_residual_gain: float
    threshold: float
    reason_code: str
    oracle_version: str
    target_state: int | None = None
    admitted: bool | None = None

    def __post_init__(self) -> None:
        if self.action not in _ALLOWED_ACTIONS:
            raise ValueError(f"action must be one of {sorted(_ALLOWED_ACTIONS)}")
        residual_gain = _finite(self.residual_gain, "residual_gain")
        normalized = _finite(self.normalized_residual_gain, "normalized_residual_gain")
        threshold = _finite(self.threshold, "threshold")
        if normalized < -1e-12:
            raise ValueError("normalized_residual_gain must be non-negative")
        if threshold < 0.0:
            raise ValueError("threshold must be non-negative")
        object.__setattr__(self, "residual_gain", residual_gain)
        object.__setattr__(self, "normalized_residual_gain", normalized)
        object.__setattr__(self, "threshold", threshold)
        object.__setattr__(self, "reason_code", _nonempty(self.reason_code, "reason_code"))
        object.__setattr__(self, "oracle_version", _nonempty(self.oracle_version, "oracle_version"))

        if self.action == "spawn":
            if self.target_state is None or int(self.target_state) != self.target_state or self.target_state < 0:
                raise ValueError("spawn decisions require a non-negative integer target_state")
            object.__setattr__(self, "target_state", int(self.target_state))
        elif self.target_state is not None:
            raise ValueError("no_spawn decisions must not carry a target_state")

        if self.admitted is not None and not isinstance(self.admitted, bool):
            raise TypeError("admitted must be a native bool or None")
        if self.action == "spawn" and self.admitted is False:
            raise ValueError("action='spawn' cannot have admitted=False")
        if self.action == "no_spawn" and self.admitted is True:
            raise ValueError("action='no_spawn' cannot have admitted=True")

    def to_payload(self) -> dict[str, Any]:
        return {
            "action": self.action,
            "target_state": self.target_state,
            "residual_gain": self.residual_gain,
            "normalized_residual_gain": self.normalized_residual_gain,
            "threshold": self.threshold,
            "reason_code": self.reason_code,
            "oracle_version": self.oracle_version,
            "admitted": self.admitted,
        }


@dataclass(frozen=True)
class SpawningRecord:
    identity: EventIdentity
    features: InvariantFeatures
    decision: OracleDecision
    candidate_kind: str
    candidate_serial: int
    schema_version: str = SCHEMA_VERSION

    def __post_init__(self) -> None:
        object.__setattr__(self, "candidate_kind", _nonempty(self.candidate_kind, "candidate_kind"))
        if int(self.candidate_serial) != self.candidate_serial or self.candidate_serial < 0:
            raise ValueError("candidate_serial must be a non-negative integer")
        object.__setattr__(self, "candidate_serial", int(self.candidate_serial))
        if self.schema_version != SCHEMA_VERSION:
            raise ValueError(f"schema_version must be {SCHEMA_VERSION}")

    def to_payload(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "identity": self.identity.to_payload(),
            "candidate": {
                "kind": self.candidate_kind,
                "serial": self.candidate_serial,
            },
            "features": self.features.to_payload(),
            "decision": self.decision.to_payload(),
        }

    @property
    def record_hash(self) -> str:
        return fingerprint(self.to_payload())


@dataclass(frozen=True)
class DatasetManifest:
    """Frozen provenance receipt for one leakage-safe dataset split."""

    split: str
    record_count: int
    record_hashes: tuple[str, ...]
    dataset_hash: str
    oracle_versions: tuple[str, ...]
    leakage_groups: tuple[str, ...]
    generation_commit: str
    frozen: bool = True
    schema_version: str = SCHEMA_VERSION

    @classmethod
    def from_records(
        cls,
        split: str,
        records: Iterable[SpawningRecord],
        generation_commit: str,
        *,
        frozen: bool = True,
    ) -> "DatasetManifest":
        if split not in _ALLOWED_SPLITS:
            raise ValueError(f"split must be one of {sorted(_ALLOWED_SPLITS)}")
        generation_commit = _nonempty(generation_commit, "generation_commit")
        materialized = tuple(records)
        hashes = tuple(sorted(record.record_hash for record in materialized))
        groups = tuple(
            sorted(
                f"{record.identity.family_id}::{record.identity.trajectory_id}"
                for record in materialized
            )
        )
        # Deduplicate groups while retaining deterministic lexical ordering.
        groups = tuple(sorted(set(groups)))
        oracle_versions = tuple(sorted({record.decision.oracle_version for record in materialized}))
        dataset_hash = fingerprint(
            {
                "schema_version": SCHEMA_VERSION,
                "split": split,
                "record_hashes": list(hashes),
                "oracle_versions": list(oracle_versions),
                "leakage_groups": list(groups),
                "generation_commit": generation_commit,
                "frozen": bool(frozen),
            }
        )
        return cls(
            split=split,
            record_count=len(materialized),
            record_hashes=hashes,
            dataset_hash=dataset_hash,
            oracle_versions=oracle_versions,
            leakage_groups=groups,
            generation_commit=generation_commit,
            frozen=bool(frozen),
        )

    def to_payload(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "split": self.split,
            "record_count": self.record_count,
            "record_hashes": list(self.record_hashes),
            "dataset_hash": self.dataset_hash,
            "oracle_versions": list(self.oracle_versions),
            "leakage_groups": list(self.leakage_groups),
            "generation_commit": self.generation_commit,
            "frozen": self.frozen,
        }


def assert_disjoint_manifests(*manifests: DatasetManifest) -> None:
    """Reject trajectory/family leakage across any pair of dataset splits."""
    seen: dict[str, str] = {}
    for manifest in manifests:
        for group in manifest.leakage_groups:
            prior = seen.get(group)
            if prior is not None and prior != manifest.split:
                raise ValueError(
                    f"leakage group {group!r} appears in both {prior!r} and {manifest.split!r}"
                )
            seen[group] = manifest.split


def hermitian_operator_invariants(operator: np.ndarray, *, atol: float = 1e-11) -> dict[str, Any]:
    """Compute invariants unchanged by unitary similarity H -> U^† H U.

    This helper is intentionally small: it provides a concrete, independently
    testable bridge from matrix-valued electronic data to admissible ML features.
    """
    matrix = np.asarray(operator, dtype=np.complex128)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError("operator must be a square matrix")
    if not np.all(np.isfinite(matrix.real)) or not np.all(np.isfinite(matrix.imag)):
        raise ValueError("operator must contain only finite entries")
    if not np.allclose(matrix, matrix.conj().T, rtol=0.0, atol=atol):
        raise ValueError("operator must be Hermitian")

    eigenvalues = np.linalg.eigvalsh(matrix).real
    gaps = np.diff(eigenvalues)
    return {
        "trace": float(np.trace(matrix).real),
        "frobenius_norm": float(np.linalg.norm(matrix, ord="fro")),
        "eigenvalues": [float(value) for value in eigenvalues],
        "adjacent_gaps": [float(value) for value in gaps],
    }


def evaluate_predictions(
    truth: Sequence[SpawningRecord],
    probabilities: Sequence[float],
    predicted_scores: Sequence[float],
    *,
    decision_threshold: float = 0.5,
) -> dict[str, float]:
    """Evaluate a future ranker without allowing it to replace the physics oracle.

    ``probabilities`` is the model's spawn probability and ``predicted_scores`` is
    its prediction of normalized residual gain.  Event-time MAE is deliberately not
    computed from flat records; trajectory-level event matching belongs in the future
    sequence evaluator and remains a required metric in the frozen contract.
    """
    if len(truth) != len(probabilities) or len(truth) != len(predicted_scores):
        raise ValueError("truth, probabilities, and predicted_scores must have equal length")
    if not truth:
        raise ValueError("evaluation requires at least one record")
    threshold = _finite(decision_threshold, "decision_threshold")
    if not 0.0 <= threshold <= 1.0:
        raise ValueError("decision_threshold must lie in [0, 1]")

    probs = np.asarray([_finite(value, "probability") for value in probabilities], dtype=float)
    if np.any((probs < 0.0) | (probs > 1.0)):
        raise ValueError("probabilities must lie in [0, 1]")
    scores = np.asarray([_finite(value, "predicted_score") for value in predicted_scores], dtype=float)
    labels = np.asarray([record.decision.action == "spawn" for record in truth], dtype=bool)
    predicted = probs >= threshold

    tp = int(np.sum(predicted & labels))
    fp = int(np.sum(predicted & ~labels))
    fn = int(np.sum(~predicted & labels))
    precision = tp / (tp + fp) if tp + fp else 1.0
    recall = tp / (tp + fn) if tp + fn else 1.0
    false_negative_rate = fn / (tp + fn) if tp + fn else 0.0
    target_scores = np.asarray(
        [record.decision.normalized_residual_gain for record in truth], dtype=float
    )

    return {
        "precision": float(precision),
        "recall": float(recall),
        "false_negative_rate": float(false_negative_rate),
        "brier_score": float(np.mean((probs - labels.astype(float)) ** 2)),
        "mean_score_error": float(np.mean(np.abs(scores - target_scores))),
    }
