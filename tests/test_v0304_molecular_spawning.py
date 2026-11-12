"""Replay and provenance controls for real Hamiltonian spawning checkpoints."""
from copy import deepcopy
from dataclasses import replace
import json

import numpy as np
import pytest

from gaussian_dynamics.curved_basis_adaptation_v292 import CurvedBasisSettingsV292
from gaussian_dynamics.molecular_midpoint_v302 import _state_payload
from gaussian_dynamics.molecular_spawn_v304 import (
    evaluate_molecular_spawn_checkpoint_v304, replay_molecular_spawn_checkpoint_v304,
    serialize_molecular_spawn_checkpoint_v304,
)
from gaussian_dynamics.pre_ml_contract_v294 import (
    V294_ORACLE_IMPLEMENTATION, V294_ORACLE_SCHEMA, SpawnOracleRecordV294, stable_hash_v294,
)
from gaussian_dynamics.v302_benchmark import free_packet_fixture_v302


@pytest.fixture(scope="module")
def model():
    return free_packet_fixture_v302()


IDENTITY = dict(family_id="free-packet", trajectory_id="trajectory-0", step=0)


@pytest.fixture(scope="module")
def checkpoint(model):
    return evaluate_molecular_spawn_checkpoint_v304(*model, **IDENTITY)


def test_every_candidate_is_recorded_without_altering_state_or_claiming_admission(model, checkpoint):
    state = model[0]
    assert checkpoint["state"] == _state_payload(state)
    assert checkpoint["candidate_count"] == len(checkpoint["candidates"]) == 4
    assert {row["candidate"]["displacement_kind"] for row in checkpoint["candidates"]} == {"position", "momentum"}
    assert {row["candidate"]["sign"] for row in checkpoint["candidates"]} == {-1, 1}
    assert checkpoint["decision_stage"] == "residual_screening"
    assert checkpoint["final_admission"] is None
    assert checkpoint["dataset_certified"] is False
    for row in checkpoint["candidates"]:
        record = row["record"]
        assert record["oracle_schema"] == V294_ORACLE_SCHEMA == "gnd-spawn-oracle-v2"
        assert record["oracle_implementation"] == V294_ORACLE_IMPLEMENTATION
        assert record["predicted_residual_reduction"] > 0
        assert record["relative_residual_reduction"] <= 1 + 1e-10
        assert record["schur_rank"] <= 2*state.nstate
        assert row["record_fingerprint"] == stable_hash_v294(record)


def test_complete_replay_is_exact(model, checkpoint):
    serialized = serialize_molecular_spawn_checkpoint_v304(checkpoint)
    assert replay_molecular_spawn_checkpoint_v304(serialized, *model[1:]) == checkpoint


def test_rejected_candidates_survive_serialization_with_original_reasons(model):
    receipt = evaluate_molecular_spawn_checkpoint_v304(
        *model, **IDENTITY, settings=CurvedBasisSettingsV292(minimum_predicted_residual_reduction=1e6),
    )
    restored = json.loads(serialize_molecular_spawn_checkpoint_v304(receipt))
    assert len(restored["candidates"]) == 4
    for row in restored["candidates"]:
        assert row["record"]["accepted_by_residual_oracle"] is False
        assert "insufficient-predicted-residual-reduction" in row["record"]["rejection_reasons"]


def test_unsealed_mutation_is_rejected(checkpoint):
    changed = deepcopy(checkpoint)
    changed["candidates"][0]["record"]["predicted_residual_reduction"] += .1
    with pytest.raises(ValueError, match="fingerprint"):
        serialize_molecular_spawn_checkpoint_v304(changed)


@pytest.mark.parametrize("corruption", ["score", "missing_candidate", "reordered_candidates", "false_admission"])
def test_rehashed_false_labels_or_incomplete_coverage_fail_recomputed_replay(model, checkpoint, corruption):
    changed = deepcopy(checkpoint)
    if corruption == "score":
        row = changed["candidates"][0]
        row["record"]["predicted_residual_reduction"] += .1
        row["record_fingerprint"] = stable_hash_v294(row["record"])
    elif corruption == "missing_candidate":
        changed["candidates"].pop()
        changed["candidate_count"] -= 1
    elif corruption == "reordered_candidates":
        changed["candidates"].reverse()
    else:
        changed["final_admission"] = True
    changed.pop("fingerprint")
    changed["fingerprint"] = stable_hash_v294(changed)
    with pytest.raises(ValueError, match="does not reproduce"):
        replay_molecular_spawn_checkpoint_v304(json.dumps(changed), *model[1:])


def test_replay_binds_operator(model, checkpoint):
    _, bridge, grid = model
    with pytest.raises(ValueError, match="Hamiltonian fingerprint"):
        replay_molecular_spawn_checkpoint_v304(
            serialize_molecular_spawn_checkpoint_v304(checkpoint), bridge,
            replace(grid, fingerprint="different-operator"),
        )


@pytest.mark.parametrize("kwargs", [
    {"family_id": ""}, {"trajectory_id": " "}, {"step": -1}, {"step": True},
    {"step": 0.5}, {"packet_ids": ()}, {"packet_ids": ("",)},
    {"tangent_step": np.nan}, {"tangent_step": 0},
])
def test_bad_identity_or_tangent_controls_fail_closed(model, kwargs):
    with pytest.raises(ValueError):
        evaluate_molecular_spawn_checkpoint_v304(*model, **{**IDENTITY, **kwargs})


def test_legacy_oracle_labels_are_not_accepted_as_current(checkpoint):
    row = deepcopy(checkpoint["candidates"][0]["record"])
    row["oracle_schema"] = "gnd-spawn-oracle-v1"
    with pytest.raises(ValueError, match="schema"):
        SpawnOracleRecordV294(**row).validate()
    row["oracle_schema"] = V294_ORACLE_SCHEMA
    row["oracle_implementation"] = "legacy-real-only"
    with pytest.raises(ValueError, match="implementation"):
        SpawnOracleRecordV294(**row).validate()
