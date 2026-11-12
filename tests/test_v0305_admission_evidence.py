"""Acceptance must be reconstructed from admission measurements and bound rows."""
from copy import deepcopy
import json
from pathlib import Path

import pytest

from gaussian_dynamics.pre_ml_contract_v294 import stable_hash_v294
from gaussian_dynamics.v305_benchmark import evaluate_admission_gates_v305


@pytest.fixture(scope="module")
def evidence():
    return json.loads((Path(__file__).resolve().parents[1] / "results/v0305_molecular_admission_evidence.json").read_text())


def seal(payload):
    payload.pop("fingerprint", None)
    payload["fingerprint"] = stable_hash_v294(payload)


def test_committed_evidence_passes_recomputed_gates(evidence):
    assert all(evaluate_admission_gates_v305(evidence).values())
    assert evidence["gates"] == evaluate_admission_gates_v305(evidence)


def test_reported_gate_flags_cannot_override_measurements(evidence):
    changed = deepcopy(evidence)
    changed["gates"] = {name: False for name in changed["gates"]}
    changed["passed"] = False
    assert all(evaluate_admission_gates_v305(changed).values())


@pytest.mark.parametrize("corruption,gate", [
    ("admission", "exactly_one_spawn"),
    ("activation", "newborn_shapes_dormant"),
    ("projection", "projection_and_energy_gates"),
    ("energy", "projection_and_energy_gates"),
    ("gauge_state", "physical_gauge_after_state"),
    ("gauge_score", "physical_gauge_scores"),
    ("gauge_features", "physical_gauge_export_features"),
    ("target", "no_invented_target_or_certification"),
    ("certification", "no_invented_target_or_certification"),
    ("export_label", "export_labels_bound_to_admission"),
    ("manifest_count", "export_manifest_counts"),
    ("split", "trajectory_splits_do_not_leak"),
    ("missing_row", "export_labels_bound_to_admission"),
    ("replay", "all_checkpoints_physically_replayed"),
])
def test_corrupted_measurements_fail_even_after_rehash(evidence, corruption, gate):
    changed = deepcopy(evidence)
    point = changed["cases"][0]
    receipt = point["fixed"]
    chosen = receipt["candidate_decisions"][receipt["selected_candidate_serial"]]
    dataset = changed["dataset"]
    row = dataset["rows"][0]
    if corruption == "admission":
        chosen["final_admitted"] = False
    elif corruption == "activation":
        receipt["after_metadata"]["active_shape_mask"][-1] = True
    elif corruption == "projection":
        chosen["attempt"]["projection"]["relative_loss"] = .5
    elif corruption == "energy":
        chosen["attempt"]["energy_jump_hartree"] = .1
    elif corruption == "gauge_state":
        point["after_wavefunction_gauge_error"] = 1e-3
    elif corruption == "gauge_score":
        point["moving"]["screening"]["candidates"][0]["record"]["predicted_residual_reduction"] += .1
    elif corruption == "gauge_features":
        row["features"]["source_coefficient_norm_squared"] += .1
    elif corruption == "target":
        row["target_state"] = 0
    elif corruption == "certification":
        dataset["training_eligible"] = True
    elif corruption == "export_label":
        row["targets"]["final_admitted"] = not row["targets"]["final_admitted"]
    elif corruption == "manifest_count":
        dataset["record_count"] += 1
    elif corruption == "split":
        row["split"] = "test" if row["split"] != "test" else "train"
    elif corruption == "missing_row":
        dataset["rows"].pop()
    else:
        dataset["exact_replay_verified"] = False
    for case in changed["cases"]:
        for frame in ("fixed", "moving"):
            seal(case[frame])
    for item in dataset["rows"]:
        seal(item)
    seal(dataset)
    assert evaluate_admission_gates_v305(changed)[gate] is False
