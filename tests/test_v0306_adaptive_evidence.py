"""Lifecycle acceptance must follow raw measurements and state bindings."""
from copy import deepcopy
import json
from pathlib import Path

import pytest

from gaussian_dynamics.pre_ml_contract_v294 import stable_hash_v294
from gaussian_dynamics.v306_benchmark import adaptive_measurements_v306, evaluate_adaptive_gates_v306


@pytest.fixture(scope="module")
def evidence():
    return json.loads((Path(__file__).resolve().parents[1] / "results/v0306_adaptive_propagation_evidence.json").read_text())


def reseal(value):
    if isinstance(value, list):
        for item in value:
            reseal(item)
    elif isinstance(value, dict):
        for item in value.values():
            reseal(item)
        if "fingerprint" in value:
            value.pop("fingerprint")
            value["fingerprint"] = stable_hash_v294(value)


def test_committed_lifecycle_evidence_passes_recomputed_gates(evidence):
    assert all(evaluate_adaptive_gates_v306(evidence).values())
    assert evidence["gates"] == evaluate_adaptive_gates_v306(evidence)
    assert evidence["measurements"] == adaptive_measurements_v306(evidence)


def test_reported_success_and_summary_fields_are_not_authoritative(evidence):
    changed = deepcopy(evidence)
    changed["passed"] = False
    changed["gates"] = {k: False for k in changed["gates"]}
    changed["measurements"] = {k: -1 for k in changed["measurements"]}
    assert all(evaluate_adaptive_gates_v306(changed).values())
    assert adaptive_measurements_v306(changed) == evidence["measurements"]


@pytest.mark.parametrize("corruption,gate", [
    ("revision", "source_revision_recorded"),
    ("schedule", "complete_explicit_schedule"),
    ("age", "state_and_metadata_continuity"),
    ("dormant_motion", "dormant_shapes_exactly_frozen"),
    ("seed", "unseeded_newborns"),
    ("growth", "both_newborns_gain_amplitude"),
    ("rank", "activation_uses_measured_gates"),
    ("speed", "activation_uses_measured_gates"),
    ("activation", "activated_shape_subsequently_moves"),
    ("admission", "two_admitted_spawns"),
    ("event_energy", "two_admitted_spawns"),
    ("dormant_topology", "topology_deferred_for_dormancy"),
    ("midpoint", "accepted_nonlinear_and_metric_residuals"),
    ("svd", "accepted_nonlinear_and_metric_residuals"),
    ("norm", "bounded_norm_and_energy_drift"),
    ("reference", "bounded_fixed_grid_reference_error"),
    ("replay", "same_runtime_physical_replay"),
    ("training", "no_training_certification"),
])
def test_corrupted_measurements_fail_after_rehash(evidence, corruption, gate):
    changed = deepcopy(evidence)
    steps = changed["trajectory"]["steps"]
    propagation = steps[0]["propagation"]
    activated = next(row for step in steps for row in step["activation"]["trials"] if row["activated"])
    if corruption == "revision":
        changed["generation_commit"] = "not-a-source-revision"
    elif corruption == "schedule":
        changed["trajectory"]["time_steps"].pop()
    elif corruption == "age":
        steps[0]["after_metadata"]["packet_ages"][0] += 3
    elif corruption == "dormant_motion":
        propagation["after"]["q"][1][0] += .1
    elif corruption == "seed":
        changed["initial_admission"]["after_state"]["coefficients_real"][-1][0] = .1
    elif corruption == "growth":
        steps[-1]["after"]["coefficients_real"][2][0] = 0
        steps[-1]["after"]["coefficients_imag"][2][0] = 0
    elif corruption == "rank":
        activated["projected_shape_metric_spectrum"][-1] = 0
    elif corruption == "speed":
        activated["trial_shape_speed"] *= 100
    elif corruption == "activation":
        activated["activated"] = False
    elif corruption == "admission":
        changed["initial_admission"]["event_kind"] = "none"
    elif corruption == "event_energy":
        row = next(r for r in changed["initial_admission"]["candidate_decisions"] if r["final_admitted"])
        row["attempt"]["energy_after_hartree"] += .1
    elif corruption == "dormant_topology":
        steps[0]["topology_status"] = "evaluated"
    elif corruption == "midpoint":
        propagation["scaled_midpoint_residual"] = 1.1
    elif corruption == "svd":
        propagation["midpoint_svd_residual"] = .1
    elif corruption == "norm":
        propagation["after_observables"]["norm"] += .01
    elif corruption == "reference":
        changed["reference"]["adaptive"]["real"][0][0] += .1
    elif corruption == "replay":
        changed["trajectory_replayed"] = False
    else:
        steps[0]["training_eligible"] = True
    reseal(changed)
    assert evaluate_adaptive_gates_v306(changed)[gate] is False


def test_unsealed_nested_receipt_is_detected(evidence):
    changed = deepcopy(evidence)
    changed["trajectory"]["steps"][0]["propagation"]["iterations"] += 1
    assert evaluate_adaptive_gates_v306(changed)["sealed_receipts"] is False


def test_baseline_error_is_computed_from_saved_wavefunction(evidence):
    changed = deepcopy(evidence)
    changed["reference"]["fixed_packet"] = deepcopy(changed["reference"]["exact"])
    assert adaptive_measurements_v306(changed)["fixed_packet_grid_exact_error"] == 0
