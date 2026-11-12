"""Ordered lifecycle steps, exact replay, and non-forgeable trajectory evidence."""
from copy import deepcopy
from dataclasses import replace
import json

import numpy as np
import pytest

from gaussian_dynamics.curved_basis_adaptation_v292 import CurvedBasisSettingsV292
from gaussian_dynamics.molecular_adaptive_v306 import (
    run_molecular_adaptive_v306, serialize_adaptive_trajectory_v306,
    replay_adaptive_trajectory_v306,
)
from gaussian_dynamics.molecular_midpoint_v302 import _state_payload
from gaussian_dynamics.pre_ml_contract_v294 import stable_hash_v294
from gaussian_dynamics.v302_benchmark import free_packet_fixture_v302

IDENTITY = dict(family_id="free", trajectory_id="adaptive")


@pytest.fixture(scope="module")
def model():
    return free_packet_fixture_v302()


@pytest.fixture(scope="module")
def trajectory(model):
    default = CurvedBasisSettingsV292()
    settings = replace(default, inherited=replace(default.inherited, maximum_packet_count=3))
    return run_molecular_adaptive_v306(*model, time_steps=[.005]*4, **IDENTITY, basis_settings=settings)


def test_propagate_age_activate_then_scheduled_topology(model, trajectory):
    assert [s.receipt["topology_status"] for s in trajectory] == ["not_scheduled", "evaluated", "not_scheduled", "deferred_dormant_shapes"]
    assert trajectory[1].receipt["topology"]["event_kind"] == "spawn"
    assert trajectory[1].receipt["after_metadata"] == {"packet_ids": ["g000000", "g000001"], "packet_ages": [2, 0], "active_shape_mask": [True, False], "next_packet_serial": 2}
    assert trajectory[-1].receipt["after_metadata"]["packet_ages"] == [4, 2]
    assert trajectory[-1].after.ngaussian == 2
    assert trajectory[-1].after.time_au == .02
    assert trajectory[0].receipt["before"] == _state_payload(model[0])
    for step in trajectory:
        r = step.receipt
        assert r["propagation"]["before"] == r["before"]
        assert r["activation"]["state_fingerprint"] == stable_hash_v294(r["propagation"]["after"])
        assert not r["dataset_certified"] and not r["training_eligible"]
        if r["topology"] is None:
            assert r["after"] == r["propagation"]["after"]


def test_full_adaptive_replay_recomputes_every_stage(model, trajectory):
    serialized = serialize_adaptive_trajectory_v306(trajectory)
    replayed = replay_adaptive_trajectory_v306(serialized, *model[1:])
    assert serialize_adaptive_trajectory_v306(replayed) == serialized


@pytest.mark.parametrize("kind", ["age", "activation", "admission"])
def test_rehashed_wrong_lifecycle_outcome_fails_replay(model, trajectory, kind):
    payload = json.loads(serialize_adaptive_trajectory_v306(trajectory))
    record = payload["steps"][0 if kind != "admission" else 1]
    if kind == "age":
        record["after_metadata"]["packet_ages"][0] += 3
    elif kind == "activation":
        record["activation"]["mask_after"][0] = False
        record["activation"].pop("fingerprint")
        record["activation"]["fingerprint"] = stable_hash_v294(record["activation"])
    else:
        record["topology"]["event_kind"] = "none"
        record["topology"].pop("fingerprint")
        record["topology"]["fingerprint"] = stable_hash_v294(record["topology"])
    record.pop("fingerprint")
    record["fingerprint"] = stable_hash_v294(record)
    payload.pop("fingerprint")
    payload["fingerprint"] = stable_hash_v294(payload)
    with pytest.raises(ValueError, match="does not reproduce"):
        replay_adaptive_trajectory_v306(json.dumps(payload), *model[1:])


def test_truncated_trace_cannot_claim_full_schedule(model, trajectory):
    payload = json.loads(serialize_adaptive_trajectory_v306(trajectory))
    payload["steps"].pop()
    payload.pop("fingerprint")
    payload["fingerprint"] = stable_hash_v294(payload)
    with pytest.raises(ValueError, match="incomplete time schedule"):
        replay_adaptive_trajectory_v306(json.dumps(payload), *model[1:])


def test_discontinuous_or_mutated_steps_cannot_be_serialized(trajectory):
    with pytest.raises(ValueError, match="contiguous"):
        serialize_adaptive_trajectory_v306([trajectory[0], trajectory[0]])
    changed = replace(trajectory[0], after=replace(trajectory[0].after, p=trajectory[0].after.p+.1))
    with pytest.raises(ValueError, match="no longer matches"):
        serialize_adaptive_trajectory_v306([changed])
    with pytest.raises(ValueError, match="empty"):
        serialize_adaptive_trajectory_v306([])


@pytest.mark.parametrize("schedule", [[], [0], [-.01], [np.nan], [np.inf], [True], [.01, 0]])
def test_invalid_forward_schedules_fail_before_evolution(model, schedule):
    before = _state_payload(model[0])
    with pytest.raises(ValueError, match="time_steps"):
        run_molecular_adaptive_v306(*model, time_steps=schedule, **IDENTITY)
    assert _state_payload(model[0]) == before


@pytest.mark.parametrize("initial_step", [-1, True, 1.2])
def test_invalid_initial_step_fails(model, initial_step):
    with pytest.raises(ValueError, match="initial_step"):
        run_molecular_adaptive_v306(*model, time_steps=[.01], initial_step=initial_step, **IDENTITY)


def test_packet_budget_control_has_no_topology_change(model):
    default = CurvedBasisSettingsV292()
    settings = replace(default, inherited=replace(default.inherited, maximum_packet_count=1, adapt_every_steps=1))
    trace = run_molecular_adaptive_v306(*model, time_steps=[.005], **IDENTITY, basis_settings=settings)
    receipt = trace[0].receipt
    assert receipt["topology"]["event_kind"] == "none"
    assert receipt["topology"]["event_reason"] == "maximum_packet_count"
    assert trace[0].after.ngaussian == 1
    assert receipt["after"] == receipt["propagation"]["after"]
