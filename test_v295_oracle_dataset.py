from dataclasses import replace

import pytest

from v294_preml_contract import EventIdentity, InvariantFeatures, OracleDecision, SpawningRecord
from v295_oracle_dataset import (
    V295_VERSION,
    audit_invariant_record_transform,
    certify_campaign,
    deterministic_group_split,
    replay_oracle_labels,
    write_jsonl,
)


def rec(family: str, trajectory: str, serial: int, *, spawn: bool = True) -> SpawningRecord:
    return SpawningRecord(
        identity=EventIdentity(
            trajectory_id=trajectory,
            family_id=family,
            event_id=f"e-{serial}",
            packet_id="p0",
            step=serial,
            time_au=0.1 * serial,
        ),
        features=InvariantFeatures(
            populations=(0.75, 0.25),
            energy_gaps=(0.02,),
            force_norms=(0.1, 0.2),
            berry_curvature_frobenius=0.03,
            gram_min_eigenvalue=0.4,
            gram_condition_number=3.0,
            residual_norm=0.5,
            novelty=0.8,
            dt_au=0.01,
        ),
        decision=OracleDecision(
            action="spawn" if spawn else "no_spawn",
            target_state=1 if spawn else None,
            residual_gain=0.2 if spawn else 0.01,
            normalized_residual_gain=0.2 if spawn else 0.01,
            threshold=0.05,
            reason_code="pass" if spawn else "below_threshold",
            oracle_version="v0.29.3-oracle",
            admitted=spawn,
        ),
        candidate_kind="principal_axis_displacement",
        candidate_serial=serial,
    )


def campaign():
    return tuple(
        rec(f"family-{i % 3}", f"traj-{i}", i, spawn=(i % 2 == 0))
        for i in range(24)
    )


def test_version():
    assert V295_VERSION == "0.29.5"


def test_split_is_deterministic_under_record_permutation():
    records = campaign()
    a = deterministic_group_split(records, salt="science")
    b = deterministic_group_split(reversed(records), salt="science")
    assert {k: [r.record_hash for r in v] for k, v in a.items()} == {
        k: [r.record_hash for r in v] for k, v in b.items()
    }


def test_split_keeps_whole_trajectory_group_together():
    records = (
        rec("f", "same", 1),
        rec("f", "same", 2, spawn=False),
        rec("f", "other", 3),
    )
    splits = deterministic_group_split(records)
    locations = {
        name for name, subset in splits.items()
        if any(r.identity.trajectory_id == "same" for r in subset)
    }
    assert len(locations) == 1
    chosen = next(iter(locations))
    assert sum(r.identity.trajectory_id == "same" for r in splits[chosen]) == 2


def test_invalid_split_fractions_fail_closed():
    with pytest.raises(ValueError, match="sum to 1"):
        deterministic_group_split(campaign(), fractions=(0.5, 0.2, 0.2, 0.2))


def test_exact_oracle_replay_passes():
    receipt = replay_oracle_labels(campaign(), lambda r: r)
    assert receipt.passed
    assert receipt.mismatch_count == 0
    assert receipt.source_hash == receipt.replay_hash


def test_oracle_replay_detects_label_change():
    def bad(r):
        if r.decision.action == "spawn":
            return replace(
                r,
                decision=replace(
                    r.decision,
                    residual_gain=r.decision.residual_gain + 1e-4,
                ),
            )
        return r

    receipt = replay_oracle_labels(campaign(), bad)
    assert not receipt.passed
    assert receipt.mismatch_count > 0
    assert receipt.source_hash != receipt.replay_hash


def test_oracle_replay_rejects_feature_mutation_even_if_label_same():
    def bad(r):
        return replace(r, features=replace(r.features, residual_norm=r.features.residual_norm + 0.01))

    receipt = replay_oracle_labels(campaign(), bad)
    assert not receipt.passed
    assert receipt.mismatch_count == len(campaign())


def test_invariance_identity_passes():
    receipt = audit_invariant_record_transform(campaign(), lambda r: r)
    assert receipt.passed


def test_invariance_audit_detects_candidate_provenance_change():
    receipt = audit_invariant_record_transform(
        campaign(), lambda r: replace(r, candidate_kind=r.candidate_kind + "_changed")
    )
    assert not receipt.passed
    assert receipt.mismatch_count == len(campaign())


def test_campaign_certification_passes_with_exact_replay_and_invariance():
    receipt = certify_campaign(
        campaign(),
        generation_commit="abc123",
        replay=lambda r: r,
        invariance_transforms={
            "electronic_gauge": lambda r: r,
            "packet_permutation": lambda r: r,
            "coordinate_transform": lambda r: r,
        },
        salt="campaign-A",
    )
    assert receipt.passed
    assert receipt.replay.passed
    assert all(x.passed for x in receipt.invariance.values())
    assert set(receipt.manifests) == {"train", "validation", "test", "audit"}


def test_campaign_hash_is_reproducible():
    kwargs = dict(
        generation_commit="abc123",
        replay=lambda r: r,
        invariance_transforms={"gauge": lambda r: r},
        salt="campaign-B",
    )
    a = certify_campaign(campaign(), **kwargs)
    b = certify_campaign(reversed(campaign()), **kwargs)
    assert a.campaign_hash == b.campaign_hash


def test_campaign_fails_when_invariance_fails():
    receipt = certify_campaign(
        campaign(),
        generation_commit="abc123",
        replay=lambda r: r,
        invariance_transforms={
            "bad_transform": lambda r: replace(
                r, features=replace(r.features, novelty=r.features.novelty - 0.01)
            )
        },
    )
    assert not receipt.passed


def test_empty_campaign_is_rejected():
    with pytest.raises(ValueError, match="at least one"):
        certify_campaign(
            (),
            generation_commit="abc123",
            replay=lambda r: r,
            invariance_transforms={},
        )


def test_jsonl_is_canonical_and_input_order_independent():
    records = campaign()[:4]
    assert write_jsonl(records) == write_jsonl(reversed(records))
    assert write_jsonl(records).endswith("\n")
