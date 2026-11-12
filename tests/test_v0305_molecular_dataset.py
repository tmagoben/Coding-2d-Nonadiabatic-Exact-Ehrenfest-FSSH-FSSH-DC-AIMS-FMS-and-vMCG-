"""Replay-verified full-spinor rows and trajectory-level split controls."""
from copy import deepcopy
from dataclasses import replace
import json

import numpy as np
import pytest

from gaussian_dynamics.molecular_admission_v305 import evaluate_molecular_admission_v305, serialize_molecular_admission_v305
from gaussian_dynamics.molecular_dataset_v305 import (
    MolecularAdmissionSourceV305, materialize_admission_dataset_v305,
    admission_dataset_jsonl_v305, trajectory_split_v305,
)
from gaussian_dynamics.pre_ml_contract_v294 import stable_hash_v294
from gaussian_dynamics.v302_benchmark import free_packet_fixture_v302

COMMIT = "a" * 40


@pytest.fixture(scope="module")
def source():
    state, bridge, grid = free_packet_fixture_v302()
    receipt = evaluate_molecular_admission_v305(state, bridge, grid, family_id="free", trajectory_id="t0", step=0).receipt
    return MolecularAdmissionSourceV305(serialize_molecular_admission_v305(receipt), bridge, grid)


@pytest.fixture(scope="module")
def dataset(source):
    return materialize_admission_dataset_v305([source], generation_commit=COMMIT)


def test_all_labels_export_without_invented_target_or_certification(dataset):
    assert dataset["record_count"] == 4 and dataset["checkpoint_count"] == 1
    assert dataset["exact_replay_verified"]
    assert not any(dataset[key] for key in ("training_eligible", "dataset_certified", "physical_invariance_certified"))
    rows = [json.loads(line) for line in admission_dataset_jsonl_v305(dataset).splitlines()]
    assert rows == dataset["rows"]
    assert sum(row["targets"]["final_admitted"] for row in rows) == 1
    for row in rows:
        assert row["target_state"] is None and row["candidate_space"] == "full_complex_electronic_spinor"
        assert row["provenance"]["generation_commit"] == COMMIT
        assert "source_coefficient_norm_squared" in row["features"]
        assert not {"tdvp_metric_condition", "tdvp_metric_rank", "predicted_residual_reduction", "source_population"}.intersection(row["features"])
        if not row["targets"]["final_admitted"]:
            assert row["targets"]["outcome"] == "not_selected"
            assert not row["targets"]["projection_attempted"]


def test_variants_share_split_and_input_order_cannot_change_dataset(source):
    variant = replace(source, variant_id="equivalent-copy")
    first = materialize_admission_dataset_v305([source, variant], generation_commit=COMMIT)
    second = materialize_admission_dataset_v305([variant, source], generation_commit=COMMIT)
    assert first == second
    assert len({row["split"] for row in first["rows"]}) == 1
    assert {row["identity"]["variant_id"] for row in first["rows"]} == {"canonical", "equivalent-copy"}


def test_duplicate_checkpoint_cannot_silently_overweight_dataset(source):
    with pytest.raises(ValueError, match="duplicate"):
        materialize_admission_dataset_v305([source, source], generation_commit=COMMIT)


def test_rehashed_false_admission_cannot_enter_dataset(source):
    receipt = json.loads(source.serialized)
    receipt["candidate_decisions"][0]["final_admitted"] = True
    receipt.pop("fingerprint")
    receipt["fingerprint"] = stable_hash_v294(receipt)
    with pytest.raises(ValueError, match="does not reproduce"):
        materialize_admission_dataset_v305([replace(source, serialized=json.dumps(receipt))], generation_commit=COMMIT)


def test_output_hash_binds_rows(dataset):
    changed = deepcopy(dataset)
    changed["rows"].pop()
    with pytest.raises(ValueError, match="fingerprint"):
        admission_dataset_jsonl_v305(changed)


@pytest.mark.parametrize("fractions", [(1, 0, 0), (1, 1, 0, 0), (np.nan, 0, 0, 0), (np.inf, 0, 0, 0), (-.1, .9, .1, .1)])
def test_invalid_splits_fail_closed(fractions):
    with pytest.raises(ValueError):
        trajectory_split_v305("family", "trajectory", fractions=fractions)


def test_split_uses_unambiguous_group_tuple():
    # These collide when joined with '::'. Across salts the tuple groups differ.
    assignments = [(trajectory_split_v305("a::b", "c", salt=str(i)), trajectory_split_v305("a", "b::c", salt=str(i))) for i in range(40)]
    assert any(left != right for left, right in assignments)
    assert trajectory_split_v305("family", "trajectory", fractions=(0, 0, 0, 1)) == "audit"


@pytest.mark.parametrize("commit", ["", "abc123", "g" * 40, "A" * 40])
def test_commit_must_be_full_provenance(source, commit):
    with pytest.raises(ValueError, match="Git SHA"):
        materialize_admission_dataset_v305([source], generation_commit=commit)


def test_empty_campaign_is_not_verified():
    with pytest.raises(ValueError, match="at least one"):
        materialize_admission_dataset_v305([], generation_commit=COMMIT)
