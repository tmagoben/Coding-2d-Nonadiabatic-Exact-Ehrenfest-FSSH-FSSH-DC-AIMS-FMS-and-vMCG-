"""Final topology decisions must follow physical and lifecycle gates."""
from copy import deepcopy
from dataclasses import replace
import json

import numpy as np
import pytest

from gaussian_dynamics import molecular_admission_v305 as admission
from gaussian_dynamics.correlated_gaussian_tdvp_v270 import CorrelatedGaussianSpinorStateV270
from gaussian_dynamics.curved_basis_adaptation_v292 import CurvedBasisSettingsV292
from gaussian_dynamics.curved_basis_adaptation_v293 import evaluate_curved_state_v293
from gaussian_dynamics.molecular_midpoint_v302 import _state_payload
from gaussian_dynamics.molecular_spawn_v304 import _checkpoint_connection
from gaussian_dynamics.pre_ml_contract_v294 import stable_hash_v294
from gaussian_dynamics.v302_benchmark import free_packet_fixture_v302

IDENTITY = dict(family_id="free", trajectory_id="admission", step=0)


@pytest.fixture(scope="module")
def model():
    return free_packet_fixture_v302()


@pytest.fixture(scope="module")
def result(model):
    return admission.evaluate_molecular_admission_v305(*model, **IDENTITY)


def settings(**values):
    default = CurvedBasisSettingsV292()
    return replace(default, inherited=replace(default.inherited, **values))


def two_packets(state, *, identical=False, tiny=False):
    return CorrelatedGaussianSpinorStateV270(
        np.vstack([state.q, state.q if identical else [[1.7]]]),
        np.vstack([state.p, state.p if identical else [[-.2]]]),
        np.repeat(state.width_matrices, 2, axis=0), np.repeat(state.chirp_matrices, 2, axis=0),
        [[1.0], [1e-8 if tiny else .5]], state.time_au,
    )


def test_spawn_is_one_event_with_dormant_newborn_and_unchanged_wavefunction(model, result):
    state, bridge, grid = model
    receipt = result.receipt
    assert receipt["event_kind"] == "spawn"
    assert receipt["screening"]["state"] == _state_payload(state)
    assert receipt["screening"]["final_admission"] is None
    assert receipt["target_state"] is None
    assert sum(row["final_admitted"] for row in receipt["candidate_decisions"]) == 1
    assert receipt["after_metadata"] == {"packet_ids": ["g000000", "g000001"], "packet_ages": [0, 0], "active_shape_mask": [True, False], "next_packet_serial": 2}
    assert result.after.time_au == state.time_au
    connection = _checkpoint_connection(bridge)
    before = evaluate_curved_state_v293(state, connection, grid.quadrature)
    after = evaluate_curved_state_v293(result.after, connection, grid.quadrature)
    assert np.linalg.norm(after - before) < 2e-12
    norm = lambda psi: np.sum(grid.quadrature.weights[:, None] * abs(psi)**2)
    energy = lambda psi: np.sum(grid.quadrature.weights[:, None] * psi.conj() * grid.apply(psi)).real / norm(psi)
    selected = receipt["candidate_decisions"][receipt["selected_candidate_serial"]]
    assert selected["attempt"]["energy_jump_hartree"] == pytest.approx(energy(after) - energy(before), abs=1e-15)
    assert not receipt["dataset_certified"] and not receipt["adaptive_trajectory_validated"]
    for row in receipt["candidate_decisions"]:
        if not row["final_admitted"]:
            assert row["outcome"] == "not_selected" and row["attempt"] is None


def test_admission_replay_recomputes_exact_after_state(model, result):
    rebuilt = admission.replay_molecular_admission_v305(admission.serialize_molecular_admission_v305(result.receipt), *model[1:])
    assert rebuilt.receipt == result.receipt
    assert _state_payload(rebuilt.after) == _state_payload(result.after)


def test_next_checkpoint_cannot_spawn_while_actual_newborn_is_dormant(model, result):
    next_result = admission.evaluate_molecular_admission_v305(
        result.after, *model[1:], **{**IDENTITY, "step": 1}, **result.receipt["after_metadata"],
    )
    assert next_result.after is result.after
    assert next_result.receipt["event_kind"] == "none"
    assert next_result.receipt["event_reason"] == "dormant_shape_pending"
    assert next_result.receipt["before_metadata"] == next_result.receipt["after_metadata"]


@pytest.mark.parametrize("kind", ["state", "admission", "energy", "activation", "missing_candidate"])
def test_rehashed_forgery_fails_physical_replay(model, result, kind):
    changed = deepcopy(result.receipt)
    selected = changed["selected_candidate_serial"]
    if kind == "state":
        changed["after_state"]["q"][0][0] += .2
    elif kind == "admission":
        changed["candidate_decisions"][selected]["final_admitted"] = False
    elif kind == "energy":
        changed["candidate_decisions"][selected]["attempt"]["energy_jump_hartree"] += .1
    elif kind == "activation":
        changed["after_metadata"]["active_shape_mask"][-1] = True
    else:
        changed["candidate_decisions"].pop()
    changed.pop("fingerprint")
    changed["fingerprint"] = stable_hash_v294(changed)
    with pytest.raises(ValueError, match="does not reproduce"):
        admission.replay_molecular_admission_v305(json.dumps(changed), *model[1:])


@pytest.mark.parametrize("overrides", [
    {"packet_ages": [-1]}, {"packet_ages": [True]}, {"packet_ages": [1.2]},
    {"packet_ids": [" padded "]}, {"packet_ids": []}, {"active_shape_mask": [1]},
    {"active_shape_mask": []}, {"next_packet_serial": 0}, {"next_packet_serial": True},
    {"next_packet_serial": -1},
])
def test_invalid_lifecycle_metadata_fails_before_physics(model, overrides):
    with pytest.raises(ValueError):
        admission.evaluate_molecular_admission_v305(*model, **IDENTITY, **overrides)


@pytest.mark.parametrize("kwargs,reason", [
    ({"active_shape_mask": [False]}, "dormant_shape_pending"),
    ({"settings": settings(maximum_packet_count=1)}, "maximum_packet_count"),
])
def test_lifecycle_blocks_are_not_residual_or_projection_rejections(model, kwargs, reason):
    result = admission.evaluate_molecular_admission_v305(*model, **IDENTITY, **kwargs)
    assert result.after is model[0]
    assert result.receipt["event_kind"] == "none"
    assert all(row["residual_accepted"] and not row["final_admitted"] and row["outcome"] == reason and row["attempt"] is None for row in result.receipt["candidate_decisions"])


def test_residual_rejections_never_attempt_projection(model):
    result = admission.evaluate_molecular_admission_v305(*model, **IDENTITY, settings=CurvedBasisSettingsV292(minimum_predicted_residual_reduction=1e6))
    assert result.after is model[0]
    assert all(row["outcome"] == "residual_rejected" and row["attempt"] is None for row in result.receipt["candidate_decisions"])


def test_real_ill_conditioned_projection_rejects_selected_candidate(model):
    result = admission.evaluate_molecular_admission_v305(*model, **IDENTITY, settings=settings(maximum_basis_condition_number=1.01))
    assert result.after is model[0]
    selected = result.receipt["candidate_decisions"][result.receipt["selected_candidate_serial"]]
    assert selected["residual_accepted"] and not selected["final_admitted"]
    assert selected["outcome"] == "projection_rejected"
    assert selected["attempt"]["projection"] is None


def test_spawn_energy_gate_rejects_injected_energy_jump(model, monkeypatch):
    energy = admission._energy
    def shifted(state, connection, grid):
        return energy(state, connection, grid) + (.01 if state.ngaussian == 2 else 0)
    monkeypatch.setattr(admission, "_energy", shifted)
    result = admission.evaluate_molecular_admission_v305(*model, **IDENTITY)
    selected = result.receipt["candidate_decisions"][result.receipt["selected_candidate_serial"]]
    assert selected["attempt"]["projection"]["relative_loss"] < 1e-20
    assert selected["outcome"] == "energy_rejected" and not selected["final_admitted"]
    assert result.after is model[0]


def test_merge_precedes_spawn_and_uses_stable_packet_ids(model):
    state, bridge, grid = model
    state = two_packets(state, identical=True)
    result = admission.evaluate_molecular_admission_v305(state, bridge, grid, **IDENTITY, packet_ids=["first", "second"])
    assert result.receipt["event_kind"] == "merge"
    assert result.receipt["removed_packet_id"] == "second"
    assert result.after.ngaussian == 1
    assert result.receipt["selected_candidate_serial"] is None
    assert all(not row["final_admitted"] for row in result.receipt["candidate_decisions"])
    assert result.receipt["removal_attempts"][0]["outcome"] == "admitted"


def test_old_tiny_packet_is_pruned_before_spawn(model):
    state, bridge, grid = model
    state = two_packets(state, tiny=True)
    result = admission.evaluate_molecular_admission_v305(state, bridge, grid, **IDENTITY, packet_ages=[100, 100])
    assert result.receipt["event_kind"] == "prune"
    assert result.receipt["removed_packet_id"] == "g000001"
    assert result.after.ngaussian == 1
    assert all(not row["final_admitted"] for row in result.receipt["candidate_decisions"])


def test_packet_permutation_preserves_removal_and_physical_state(model):
    state, bridge, grid = model
    state = two_packets(state, identical=True)
    a = admission.evaluate_molecular_admission_v305(state, bridge, grid, **IDENTITY, packet_ids=["first", "second"])
    permuted = replace(state, q=state.q[::-1], p=state.p[::-1], width_matrices=state.width_matrices[::-1], chirp_matrices=state.chirp_matrices[::-1], coefficients=state.coefficients[::-1])
    b = admission.evaluate_molecular_admission_v305(permuted, bridge, grid, **IDENTITY, packet_ids=["second", "first"])
    assert a.receipt["removed_packet_id"] == b.receipt["removed_packet_id"] == "second"
    assert a.receipt["after_metadata"] == b.receipt["after_metadata"]
    assert np.allclose(a.after.coefficients, b.after.coefficients, atol=2e-12, rtol=0)


def test_pruning_energy_jump_is_measured_and_rejected(model):
    state, bridge, grid = model
    state = two_packets(state)
    controls = settings(prune_population_threshold=.3, maximum_projection_loss=1.0,
                        maximum_event_energy_jump_hartree=1e-12, maximum_packet_count=2)
    result = admission.evaluate_molecular_admission_v305(state, bridge, grid, **IDENTITY,
                                                       packet_ages=[100, 100], settings=controls)
    assert result.after is state
    assert result.receipt["event_kind"] == "none"
    attempt = result.receipt["removal_attempts"][0]
    assert attempt["kind"] == "prune" and attempt["outcome"] == "energy_rejected"
    assert attempt["projection"]["relative_loss"] < 1.0
    assert abs(attempt["energy_jump_hartree"]) > 1e-6
