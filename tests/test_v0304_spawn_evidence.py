"""Numerical acceptance controls; CI separately regenerates the SOC snapshots."""
import json
from pathlib import Path

import numpy as np
import pytest

from gaussian_dynamics.pre_ml_contract_v294 import stable_hash_v294
from gaussian_dynamics.v304_benchmark import SCHEMA_V304_EVIDENCE, evaluate_spawn_gates_v304


@pytest.fixture
def evidence():
    return json.loads((Path(__file__).resolve().parents[1] / "results/v0304_molecular_spawning_evidence.json").read_text())


def test_committed_campaign_is_sealed_and_scope_is_explicit(evidence):
    digest = evidence.pop("fingerprint")
    assert digest == stable_hash_v294(evidence)
    assert evidence["schema"] == SCHEMA_V304_EVIDENCE
    assert evidence["gates"] == evaluate_spawn_gates_v304(evidence)
    assert evidence["passed"] is True
    assert all(evidence["gates"].values())
    assert evidence["dataset_certified"] is False
    assert evidence["live_molecular_soc_validated"] is False
    assert evidence["propagation_receipt"]["scaled_midpoint_residual"] <= 1


@pytest.mark.parametrize("value", [np.nan, np.inf, 1.0])
@pytest.mark.parametrize("field,gate", [
    ("moving_gain", "physical_gauge_scores"),
    ("moving_relative_gain", "physical_gauge_scores"),
    ("fixed_reference_gain", "independent_complex_lstsq"),
])
def test_bad_scores_cannot_be_hidden_by_passing_flags(evidence, value, field, gate):
    evidence["snapshots"][0]["comparisons"][0][field] = value
    evidence["gates"] = {name: True for name in evidence["gates"]}
    assert not evaluate_spawn_gates_v304(evidence)[gate]


def test_discrete_gauge_changes_and_incomplete_candidate_coverage_fail(evidence):
    row = evidence["snapshots"][0]["comparisons"][0]
    row["moving_accepted"] = not row["fixed_accepted"]
    row["moving_rank"] += 1
    evidence["snapshots"][0]["moving"]["candidates"].pop()
    gates = evaluate_spawn_gates_v304(evidence)
    assert not gates["physical_gauge_decisions"]
    assert not gates["physical_gauge_ranks"]
    assert not gates["complete_candidate_coverage"]


def test_lost_rejection_reason_failed_replay_and_false_admission_fail(evidence):
    evidence["rejected_checkpoint"]["candidates"][0]["record"]["rejection_reasons"] = []
    evidence["same_runtime_exact_replay"] = False
    evidence["snapshots"][0]["fixed"]["final_admission"] = True
    gates = evaluate_spawn_gates_v304(evidence)
    assert not gates["residual_rejections_recorded"]
    assert not gates["same_runtime_exact_replay"]
    assert not gates["residual_labels_are_not_final_admission"]


def test_electronic_only_action_cannot_pass_as_full_hamiltonian(evidence):
    evidence["snapshots"][0]["kinetic_action_norm"] = 0.0
    assert not evaluate_spawn_gates_v304(evidence)["full_kinetic_action"]


def test_reported_comparisons_must_match_the_underlying_candidate_records(evidence):
    # Even self-consistent comparison/reference columns cannot silently replace
    # the independently serialized physical checkpoint.
    comparison = evidence["snapshots"][0]["comparisons"][0]
    for name in ("fixed_gain", "moving_gain", "fixed_reference_gain", "moving_reference_gain"):
        comparison[name] = .25
    assert not evaluate_spawn_gates_v304(evidence)["comparison_records_match_checkpoints"]


@pytest.mark.parametrize("value", [-.1, 1.1])
def test_gain_fraction_outside_physical_bounds_fails(evidence, value):
    evidence["snapshots"][0]["comparisons"][0]["fixed_relative_gain"] = value
    assert not evaluate_spawn_gates_v304(evidence)["bounded_residual_gains"]
