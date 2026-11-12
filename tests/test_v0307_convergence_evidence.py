"""Numerical acceptance is recomputed from physical fields and step diagnostics."""
from copy import deepcopy
import json
from pathlib import Path

import numpy as np
import pytest

from gaussian_dynamics.pre_ml_contract_v294 import stable_hash_v294
from gaussian_dynamics.v307_benchmark import convergence_measurements_v307, evaluate_convergence_gates_v307, run_v307_convergence_campaign


@pytest.fixture(scope="module")
def evidence():
    return json.loads((Path(__file__).resolve().parents[1]/"results/v0307_adaptive_soc_convergence_evidence.json").read_text())


def reseal(value):
    if isinstance(value,list):
        for x in value:
            reseal(x)
    elif isinstance(value,dict):
        for x in value.values():
            reseal(x)
        if "fingerprint" in value:
            value.pop("fingerprint")
            value["fingerprint"] = stable_hash_v294(value)


def test_committed_campaign_passes_all_recomputed_gates(evidence):
    assert all(evaluate_convergence_gates_v307(evidence).values())
    assert evidence["gates"] == evaluate_convergence_gates_v307(evidence)
    assert evidence["measurements"] == convergence_measurements_v307(evidence)
    payload = dict(evidence)
    digest = payload.pop("fingerprint")
    assert stable_hash_v294(payload) == digest


def test_saved_success_and_summary_fields_are_not_authoritative(evidence):
    changed = deepcopy(evidence)
    changed["passed"] = False
    changed["gates"] = {}
    changed["measurements"] = {}
    assert all(evaluate_convergence_gates_v307(changed).values())
    assert convergence_measurements_v307(changed) == evidence["measurements"]


@pytest.mark.parametrize("kind,gate",[
    ("revision","source_revision_recorded"),
    ("coverage","complete_two_case_two_gauge_campaign"),
    ("schedule","fixed_physical_checkpoint_schedule"),
    ("age","state_clock_and_metadata_continuity"),
    ("microstep","state_clock_and_metadata_continuity"),
    ("seed","unseeded_admitted_births"),
    ("energy_jump","unseeded_admitted_births"),
    ("dormant","dormant_shapes_exactly_frozen"),
    ("activation_speed","measured_activation_gates"),
    ("activation_rank","measured_activation_gates"),
    ("midpoint","accepted_midpoint_metric_and_null_residuals"),
    ("svd","accepted_midpoint_metric_and_null_residuals"),
    ("null","accepted_midpoint_metric_and_null_residuals"),
    ("drift","bounded_norm_and_energy_drift"),
    ("order","second_order_fixed_checkpoint_convergence"),
    ("gauge","physical_gauge_error_converges"),
    ("reference","bounded_finite_grid_reference_error"),
    ("off","physical_soc_transfer_and_disabled_control"),
    ("cadence","checkpoint_cadence_has_controlled_microsteps"),
    ("bracket","activation_bracket_narrows_with_checkpoint_cadence"),
    ("replay","selected_complete_physical_replays"),
    ("certification","no_training_or_continuum_certification"),
])
def test_rehashed_corrupted_measurements_fail(evidence,kind,gate):
    changed = deepcopy(evidence)
    frame = changed["cases"][0]["frames"]["fixed"]
    run = frame["levels"][0]
    row = run["checkpoints"][0]
    propagation = row["propagations"][0]
    activation = next(x for c in run["checkpoints"] for x in c["activation"]["trials"] if x["activated"])
    if kind == "revision":
        changed["generation_commit"] = "unknown"
    elif kind == "coverage":
        changed["cases"][1]["q0"] = 0
    elif kind == "schedule":
        run["checkpoint_times"][0] += .001
    elif kind == "age":
        row["after_metadata"]["packet_ages"][0] += 1
    elif kind == "microstep":
        propagation["after_time"] += .001
    elif kind == "seed":
        frame["initial_admission"]["after_state"]["coefficients_real"][1][0] = .01
    elif kind == "energy_jump":
        selected = next(x["attempt"] for x in frame["initial_admission"]["candidate_decisions"] if x["final_admitted"])
        selected["energy_after_hartree"] += .01
    elif kind == "dormant":
        row["after"]["q"][1][0] += .01
    elif kind == "activation_speed":
        activation["trial_shape_speed"] *= 100
    elif kind == "activation_rank":
        activation["projected_shape_metric_spectrum"][-1] = 0
    elif kind == "midpoint":
        propagation["scaled_midpoint_residual"] = 1.01
    elif kind == "svd":
        propagation["midpoint_svd_residual"] = .1
    elif kind == "null":
        propagation["midpoint_null_rhs_relative"] = .1
    elif kind == "drift":
        propagation["after_observables"]["norm"] += .001
    elif kind == "order":
        run["final_field"] = deepcopy(frame["levels"][1]["final_field"])
    elif kind == "gauge":
        changed["cases"][0]["frames"]["moving"]["levels"][2]["final_field"]["real"][0][0] += .001
    elif kind == "reference":
        run["exact_field"]["real"][0][0] += 1
    elif kind == "off":
        changed["soc_off"]["final_field"]["real"][0][1] = .1
    elif kind == "cadence":
        changed["cadence"][-1]["checkpoints"][0]["propagations"][0]["dt"] *= 2
    elif kind == "bracket":
        for c in changed["cadence"][-1]["checkpoints"]:
            for x in c["activation"]["trials"]:
                x["activated"] = False
    elif kind == "replay":
        run["full_replay_matches"] = False
    else:
        changed["continuum_certified"] = True
    reseal(changed)
    assert evaluate_convergence_gates_v307(changed)[gate] is False


def test_checkpoint_brackets_are_sampling_intervals_not_event_root_claims(evidence):
    measured = convergence_measurements_v307(evidence)
    np.testing.assert_allclose([x["activation_times"][0] for x in measured["cadence"]["levels"]],[.04,.04,.035,.0325],rtol=0,atol=1e-14)
    assert "observed_order" not in measured["cadence"]
    # The unchanged first event time produces a plateau, not a fabricated rate.
    assert measured["cadence"]["successive_endpoint_differences"][0] < 1e-7


def test_unsealed_activation_is_detected(evidence):
    changed = deepcopy(evidence)
    changed["cases"][0]["frames"]["fixed"]["levels"][0]["checkpoints"][0]["activation"]["tangent_step"] *= 2
    assert not evaluate_convergence_gates_v307(changed)["sealed_admissions_and_activations"]


@pytest.mark.parametrize("revision",[None,"main","abc",True])
def test_missing_source_revision_fails_before_any_computation(revision):
    with pytest.raises(ValueError,match="source commit SHA"):
        run_v307_convergence_campaign(generation_commit=revision)
