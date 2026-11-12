"""Physical adaptation clocks remain independent of the integration mesh."""
from dataclasses import replace
import json

import numpy as np
import pytest

from gaussian_dynamics.analytic_soc_bridge_v307 import analytic_soc_fixture_v307
from gaussian_dynamics.curved_basis_adaptation_v292 import CurvedBasisSettingsV292
from gaussian_dynamics.molecular_admission_v305 import evaluate_molecular_admission_v305
from gaussian_dynamics.molecular_checkpoints_v307 import (
    run_molecular_checkpoints_v307, serialize_checkpoint_trajectory_v307,
    replay_checkpoint_trajectory_v307,
)
from gaussian_dynamics.molecular_midpoint_v302 import _state_payload
from gaussian_dynamics.pre_ml_contract_v294 import stable_hash_v294

IDENTITY = dict(family_id="soc", trajectory_id="physical-checkpoints")


@pytest.fixture(scope="module")
def model():
    state, bridge, grid, _, _ = analytic_soc_fixture_v307()
    default = CurvedBasisSettingsV292()
    settings = replace(default, inherited=replace(default.inherited, maximum_packet_count=2, adapt_every_steps=1))
    initial = evaluate_molecular_admission_v305(state, bridge, grid, **IDENTITY, step=0, settings=settings)
    return initial, bridge, grid, settings


def run(model, **kwargs):
    initial, bridge, grid, settings = model
    return run_molecular_checkpoints_v307(initial.after, bridge, grid, **IDENTITY,
                                         **initial.receipt["after_metadata"], basis_settings=settings, **kwargs)


@pytest.fixture(scope="module")
def trajectory(model):
    return run(model, checkpoint_times=[.02, .04], substeps=2)


def test_mesh_refinement_preserves_checkpoint_age_and_activation_clock(model, trajectory):
    coarse = run(model, checkpoint_times=[.02, .04], substeps=1)
    assert model[0].receipt["event_kind"] == "spawn"
    assert np.linalg.norm(model[0].after.coefficients[-1])**2 < 1e-28  # projection roundoff, no amplitude seed
    for a, b in zip(coarse, trajectory):
        assert a.after.time_au == b.after.time_au == a.receipt["checkpoint_time"]
        assert a.receipt["after_metadata"] == b.receipt["after_metadata"]
        assert len(a.receipt["propagations"]) == 1
        assert len(b.receipt["propagations"]) == 2
    assert trajectory[0].receipt["topology_status"] == "deferred_dormant_shapes"
    assert trajectory[0].receipt["activation"]["trials"][0]["reason"] == "shape_velocity_amplification"
    assert trajectory[1].receipt["after_metadata"]["active_shape_mask"] == [True, True]
    assert trajectory[1].receipt["after_metadata"]["packet_ages"] == [2, 2]
    assert trajectory[1].receipt["topology"]["event_reason"] == "maximum_packet_count"


def test_every_microstep_is_contiguous_and_dormant_shapes_are_exactly_fixed(trajectory):
    for result in trajectory:
        r = result.receipt
        state = r["before"]
        for propagation in r["propagations"]:
            assert propagation["before"] == state
            for name in ("q", "p", "width_matrices", "chirp_matrices"):
                assert propagation["before"][name][1] == propagation["after"][name][1]
            assert propagation["scaled_midpoint_residual"] <= 1
            state = propagation["after"]
        assert r["activation"]["state_fingerprint"] == stable_hash_v294(state)
        assert not r["dataset_certified"] and not r["training_eligible"]


def test_full_checkpoint_replay_recomputes_microsteps_and_lifecycle(model, trajectory):
    serialized = serialize_checkpoint_trajectory_v307(trajectory)
    assert serialize_checkpoint_trajectory_v307(replay_checkpoint_trajectory_v307(serialized, *model[1:3])) == serialized


@pytest.mark.parametrize("kind", ["age", "substeps", "activation", "residual"])
def test_rehashed_wrong_checkpoint_fails_physical_replay(model, trajectory, kind):
    payload = json.loads(serialize_checkpoint_trajectory_v307(trajectory))
    row = payload["steps"][0]
    if kind == "age":
        row["after_metadata"]["packet_ages"][0] += 1
    elif kind == "substeps":
        row["substeps"] += 1
    elif kind == "activation":
        row["activation"]["mask_after"][1] = True
    else:
        row["propagations"][0]["scaled_midpoint_residual"] = 0
    row.pop("fingerprint")
    row["fingerprint"] = stable_hash_v294(row)
    payload.pop("fingerprint")
    payload["fingerprint"] = stable_hash_v294(payload)
    with pytest.raises(ValueError, match="does not reproduce"):
        replay_checkpoint_trajectory_v307(json.dumps(payload), *model[1:3])


def test_truncated_schedule_and_discontinuous_state_cannot_serialize_or_replay(model, trajectory):
    payload = json.loads(serialize_checkpoint_trajectory_v307(trajectory))
    payload["steps"].pop()
    payload.pop("fingerprint")
    payload["fingerprint"] = stable_hash_v294(payload)
    with pytest.raises(ValueError, match="incomplete checkpoint time schedule"):
        replay_checkpoint_trajectory_v307(json.dumps(payload), *model[1:3])
    with pytest.raises(ValueError, match="contiguous"):
        serialize_checkpoint_trajectory_v307([trajectory[0], trajectory[0]])
    with pytest.raises(ValueError, match="no longer match"):
        serialize_checkpoint_trajectory_v307([replace(trajectory[0], after=replace(trajectory[0].after, p=trajectory[0].after.p+.1))])
    with pytest.raises(ValueError, match="empty"):
        serialize_checkpoint_trajectory_v307([])


@pytest.mark.parametrize("times", [[], [0], [-.01], [np.nan], [np.inf], [True], [.02, .01]])
def test_invalid_schedule_fails_before_evolution(model, times):
    before = _state_payload(model[0].after)
    with pytest.raises(ValueError, match="checkpoint_times"):
        run(model, checkpoint_times=times, substeps=1)
    assert before == _state_payload(model[0].after)


@pytest.mark.parametrize("count", [0, -1, True, 1.5])
def test_invalid_mesh_fails_before_evolution(model, count):
    with pytest.raises(ValueError, match="substeps"):
        run(model, checkpoint_times=[.02], substeps=count)


def test_failed_microstep_does_not_advance_input_or_publish_checkpoint(model, monkeypatch):
    import gaussian_dynamics.molecular_checkpoints_v307 as module
    original = module.molecular_active_midpoint_step_v306
    calls, published = [], []
    before = _state_payload(model[0].after)
    def fail_second(*args, **kwargs):
        calls.append(1)
        if len(calls) == 2:
            raise ValueError("injected solver failure")
        return original(*args, **kwargs)
    monkeypatch.setattr(module, "molecular_active_midpoint_step_v306", fail_second)
    with pytest.raises(ValueError, match="injected solver failure"):
        run(model, checkpoint_times=[.02], substeps=2, progress=published.append)
    assert calls == [1, 1] and not published
    assert before == _state_payload(model[0].after)


def test_checkpoint_cadence_is_a_separate_physical_policy(model, trajectory):
    finer = run(model, checkpoint_times=[.01, .02, .03, .04], substeps=1)
    assert finer[-1].receipt["after_metadata"]["packet_ages"] == [4, 4]
    assert trajectory[-1].receipt["after_metadata"]["packet_ages"] == [2, 2]
    assert finer[-1].receipt["after_metadata"]["active_shape_mask"] == [True, True]
