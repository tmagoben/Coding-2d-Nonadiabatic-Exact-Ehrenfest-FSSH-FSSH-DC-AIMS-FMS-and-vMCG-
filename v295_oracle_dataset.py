"""v0.29.5 deterministic oracle-replay and dataset-certification layer.

This module connects the frozen v0.29.4 record contract to trajectory campaigns.
It deliberately accepts already-evaluated physics-oracle records: no learned model
can create or alter an OracleDecision here.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Iterable, Mapping, Sequence

from v294_preml_contract import (
    DatasetManifest,
    SpawningRecord,
    assert_disjoint_manifests,
    canonical_json,
    fingerprint,
)

V295_VERSION = "0.29.5"
_SPLITS = ("train", "validation", "test", "audit")


def _record_key(record: SpawningRecord) -> tuple[str, str, int, int]:
    i = record.identity
    return (i.family_id, i.trajectory_id, i.step, record.candidate_serial)


def _group_key(record: SpawningRecord) -> str:
    i = record.identity
    return f"{i.family_id}::{i.trajectory_id}"


def deterministic_group_split(
    records: Iterable[SpawningRecord],
    *,
    salt: str = "v0.29.5",
    fractions: Sequence[float] = (0.70, 0.15, 0.10, 0.05),
) -> dict[str, tuple[SpawningRecord, ...]]:
    """Assign complete family/trajectory groups to deterministic dataset splits."""
    records = tuple(records)
    if len(fractions) != 4 or any(float(x) < 0.0 for x in fractions):
        raise ValueError("fractions must contain four non-negative values")
    total = float(sum(fractions))
    if abs(total - 1.0) > 1e-12:
        raise ValueError("fractions must sum to 1")
    salt = str(salt).strip()
    if not salt:
        raise ValueError("salt must be non-empty")

    cuts = []
    running = 0.0
    for value in fractions[:-1]:
        running += float(value)
        cuts.append(running)

    grouped: dict[str, list[SpawningRecord]] = {}
    for record in records:
        grouped.setdefault(_group_key(record), []).append(record)

    out: dict[str, list[SpawningRecord]] = {name: [] for name in _SPLITS}
    for group, members in sorted(grouped.items()):
        # 64 bits are ample and keep the mapping exactly reproducible.
        u = int(fingerprint({"salt": salt, "group": group})[:16], 16) / float(16**16)
        idx = 0
        while idx < len(cuts) and u >= cuts[idx]:
            idx += 1
        out[_SPLITS[idx]].extend(sorted(members, key=_record_key))
    return {name: tuple(values) for name, values in out.items()}


@dataclass(frozen=True)
class ReplayReceipt:
    record_count: int
    checked_count: int
    mismatch_count: int
    source_hash: str
    replay_hash: str
    passed: bool

    def to_payload(self) -> dict[str, object]:
        return {
            "version": V295_VERSION,
            "record_count": self.record_count,
            "checked_count": self.checked_count,
            "mismatch_count": self.mismatch_count,
            "source_hash": self.source_hash,
            "replay_hash": self.replay_hash,
            "passed": self.passed,
        }


def replay_oracle_labels(
    records: Iterable[SpawningRecord],
    replay: Callable[[SpawningRecord], SpawningRecord],
) -> ReplayReceipt:
    """Re-evaluate records and require exact decision-payload reproduction.

    The callback must return a complete SpawningRecord with the same immutable
    identity, candidate provenance, and invariant features. This makes accidental
    feature mutation or relabeling visible rather than checking action alone.
    """
    source = tuple(sorted(records, key=_record_key))
    replayed: list[SpawningRecord] = []
    mismatches = 0
    for record in source:
        candidate = replay(record)
        if not isinstance(candidate, SpawningRecord):
            raise TypeError("replay callback must return SpawningRecord")
        replayed.append(candidate)
        if (
            candidate.identity != record.identity
            or candidate.candidate_kind != record.candidate_kind
            or candidate.candidate_serial != record.candidate_serial
            or candidate.features != record.features
            or candidate.decision.to_payload() != record.decision.to_payload()
        ):
            mismatches += 1

    source_hash = fingerprint([record.to_payload() for record in source])
    replay_hash = fingerprint([record.to_payload() for record in replayed])
    return ReplayReceipt(
        record_count=len(source),
        checked_count=len(replayed),
        mismatch_count=mismatches,
        source_hash=source_hash,
        replay_hash=replay_hash,
        passed=(mismatches == 0 and source_hash == replay_hash),
    )


@dataclass(frozen=True)
class InvarianceReceipt:
    record_count: int
    transformed_count: int
    mismatch_count: int
    reference_hash: str
    transformed_hash: str
    passed: bool

    def to_payload(self) -> dict[str, object]:
        return {
            "version": V295_VERSION,
            "record_count": self.record_count,
            "transformed_count": self.transformed_count,
            "mismatch_count": self.mismatch_count,
            "reference_hash": self.reference_hash,
            "transformed_hash": self.transformed_hash,
            "passed": self.passed,
        }


def audit_invariant_record_transform(
    records: Iterable[SpawningRecord],
    transform: Callable[[SpawningRecord], SpawningRecord],
) -> InvarianceReceipt:
    """Require a physical symmetry transform to leave the ML-facing record unchanged."""
    reference = tuple(sorted(records, key=_record_key))
    transformed: list[SpawningRecord] = []
    mismatches = 0
    for record in reference:
        candidate = transform(record)
        if not isinstance(candidate, SpawningRecord):
            raise TypeError("transform callback must return SpawningRecord")
        transformed.append(candidate)
        if candidate.to_payload() != record.to_payload():
            mismatches += 1
    reference_hash = fingerprint([record.to_payload() for record in reference])
    transformed_hash = fingerprint([record.to_payload() for record in transformed])
    return InvarianceReceipt(
        record_count=len(reference),
        transformed_count=len(transformed),
        mismatch_count=mismatches,
        reference_hash=reference_hash,
        transformed_hash=transformed_hash,
        passed=(mismatches == 0 and reference_hash == transformed_hash),
    )


@dataclass(frozen=True)
class CampaignReceipt:
    generation_commit: str
    manifests: Mapping[str, DatasetManifest]
    replay: ReplayReceipt
    invariance: Mapping[str, InvarianceReceipt]
    campaign_hash: str
    passed: bool

    def to_payload(self) -> dict[str, object]:
        return {
            "version": V295_VERSION,
            "generation_commit": self.generation_commit,
            "manifests": {k: v.to_payload() for k, v in sorted(self.manifests.items())},
            "replay": self.replay.to_payload(),
            "invariance": {k: v.to_payload() for k, v in sorted(self.invariance.items())},
            "campaign_hash": self.campaign_hash,
            "passed": self.passed,
        }


def certify_campaign(
    records: Iterable[SpawningRecord],
    *,
    generation_commit: str,
    replay: Callable[[SpawningRecord], SpawningRecord],
    invariance_transforms: Mapping[str, Callable[[SpawningRecord], SpawningRecord]],
    salt: str = "v0.29.5",
    fractions: Sequence[float] = (0.70, 0.15, 0.10, 0.05),
) -> CampaignReceipt:
    """Build manifests and fail closed unless replay and invariance audits pass."""
    generation_commit = str(generation_commit).strip()
    if not generation_commit:
        raise ValueError("generation_commit must be non-empty")
    records = tuple(records)
    if not records:
        raise ValueError("campaign requires at least one record")

    splits = deterministic_group_split(records, salt=salt, fractions=fractions)
    manifests = {
        name: DatasetManifest.from_records(name, subset, generation_commit)
        for name, subset in splits.items()
    }
    assert_disjoint_manifests(*manifests.values())

    replay_receipt = replay_oracle_labels(records, replay)
    invariance = {
        name: audit_invariant_record_transform(records, transform)
        for name, transform in sorted(invariance_transforms.items())
    }
    passed = replay_receipt.passed and all(receipt.passed for receipt in invariance.values())

    payload = {
        "version": V295_VERSION,
        "generation_commit": generation_commit,
        "manifests": {k: v.to_payload() for k, v in sorted(manifests.items())},
        "replay": replay_receipt.to_payload(),
        "invariance": {k: v.to_payload() for k, v in sorted(invariance.items())},
        "passed": passed,
    }
    return CampaignReceipt(
        generation_commit=generation_commit,
        manifests=manifests,
        replay=replay_receipt,
        invariance=invariance,
        campaign_hash=fingerprint(payload),
        passed=passed,
    )


def write_jsonl(records: Iterable[SpawningRecord]) -> str:
    """Canonical JSONL serialization ordered independently of input iteration."""
    ordered = sorted(records, key=_record_key)
    return "".join(canonical_json(record.to_payload()) + "\n" for record in ordered)
