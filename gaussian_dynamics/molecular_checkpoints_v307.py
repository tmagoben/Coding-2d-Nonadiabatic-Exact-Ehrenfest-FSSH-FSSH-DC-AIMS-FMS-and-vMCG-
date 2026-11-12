"""Adaptive topology at fixed physical times, independently of integration mesh."""
from dataclasses import asdict, dataclass
import json

import numpy as np

from .curved_basis_adaptation_v292 import CurvedBasisSettingsV292
from .molecular_active_v306 import molecular_active_midpoint_step_v306, activate_molecular_shapes_v306
from .molecular_admission_v305 import _metadata, evaluate_molecular_admission_v305
from .molecular_midpoint_v302 import MolecularMidpointSettingsV302, _state_payload, _state_from_payload
from .molecular_spawn_v304 import _name, _restore_settings
from .pre_ml_contract_v294 import stable_hash_v294

SCHEMA_V307_CHECKPOINT = "gnd-physical-time-adaptive-checkpoints-v0.30.7-dev"
POLICY_V307 = "integrate to checkpoint; age once; activate; scheduled topology only with all shapes active"


@dataclass(frozen=True)
class MolecularCheckpointV307:
    before: object
    after: object
    receipt: dict


def molecular_checkpoint_step_v307(
    state, bridge, hamiltonian, end_time, *, substeps, family_id, trajectory_id, step,
    packet_ids=None, packet_ages=None, active_shape_mask=None, next_packet_serial=None,
    basis_settings=CurvedBasisSettingsV292(), midpoint_settings=MolecularMidpointSettingsV302(),
):
    """Integrate a fixed mask, then execute exactly one physical-time checkpoint.

    Ages count checkpoints, never internal integration substeps. Failed substeps
    return no accepted checkpoint and leave the caller's input state untouched.
    """
    state.validate(require_normalized=False)
    basis_settings.validate()
    midpoint_settings.validate()
    if isinstance(end_time, (bool, np.bool_)) or not np.isfinite(end_time) or end_time <= state.time_au:
        raise ValueError("end_time must be finite and later than the state time")
    if type(substeps) is not int or substeps < 1:
        raise ValueError("substeps must be a positive integer")
    if type(step) is not int or step < 1:
        raise ValueError("checkpoint step must be a positive integer")
    identity = dict(family_id=_name(family_id, "family_id"), trajectory_id=_name(trajectory_id, "trajectory_id"), step=step)
    ids, ages, mask, serial = _metadata(state, packet_ids, packet_ages, active_shape_mask, next_packet_serial)
    before_metadata = dict(packet_ids=list(ids), packet_ages=list(ages), active_shape_mask=list(mask), next_packet_serial=serial)
    times = np.linspace(state.time_au, float(end_time), substeps+1)
    if np.any(np.diff(times) <= 0):
        raise ValueError("checkpoint substeps are not representable")
    current, propagations = state, []
    for target in times[1:]:
        propagated = molecular_active_midpoint_step_v306(
            current, bridge, hamiltonian, float(target-current.time_au), active_shape_mask=mask, settings=midpoint_settings,
        )
        current = propagated.after
        propagations.append(propagated.receipt)
    active, activation = activate_molecular_shapes_v306(
        current, bridge, hamiltonian, mask, packet_ids=ids, settings=basis_settings, tangent_step=midpoint_settings.tangent_step,
    )
    after_metadata = dict(packet_ids=list(ids), packet_ages=[age+1 for age in ages], active_shape_mask=active.tolist(), next_packet_serial=serial)
    topology = None
    if step % basis_settings.inherited.adapt_every_steps:
        status = "not_scheduled"
    elif not np.all(active):
        status = "deferred_dormant_shapes"
    else:
        status = "evaluated"
        admission = evaluate_molecular_admission_v305(
            current, bridge, hamiltonian, **identity, **after_metadata,
            settings=basis_settings, tangent_step=midpoint_settings.tangent_step,
        )
        current, topology = admission.after, admission.receipt
        after_metadata = topology["after_metadata"]
    receipt = {
        "schema": SCHEMA_V307_CHECKPOINT, "policy": POLICY_V307, "identity": identity,
        "checkpoint_time": float(end_time), "substeps": substeps,
        "basis_settings": asdict(basis_settings), "midpoint_settings": asdict(midpoint_settings),
        "before": _state_payload(state), "before_metadata": before_metadata,
        "propagations": propagations, "activation": activation,
        "topology_status": status, "topology": topology,
        "after": _state_payload(current), "after_metadata": after_metadata,
        "dataset_certified": False, "training_eligible": False,
    }
    receipt["fingerprint"] = stable_hash_v294(receipt)
    return MolecularCheckpointV307(state, current, receipt)


def run_molecular_checkpoints_v307(
    state, bridge, hamiltonian, *, checkpoint_times, substeps, family_id, trajectory_id, initial_step=0,
    packet_ids=None, packet_ages=None, active_shape_mask=None, next_packet_serial=None,
    basis_settings=CurvedBasisSettingsV292(), midpoint_settings=MolecularMidpointSettingsV302(), progress=None,
):
    times = tuple(checkpoint_times)
    if not times or any(isinstance(t, (bool, np.bool_)) or not np.isfinite(t) for t in times):
        raise ValueError("checkpoint_times must be a nonempty finite real schedule")
    if np.any(np.diff([state.time_au, *times]) <= 0):
        raise ValueError("checkpoint_times must be strictly increasing after the state time")
    if type(substeps) is not int or substeps < 1:
        raise ValueError("substeps must be a positive integer")
    if type(initial_step) is not int or initial_step < 0:
        raise ValueError("initial_step must be a nonnegative integer")
    metadata = dict(packet_ids=packet_ids, packet_ages=packet_ages, active_shape_mask=active_shape_mask, next_packet_serial=next_packet_serial)
    trajectory = []
    for index, end_time in enumerate(times, initial_step+1):
        result = molecular_checkpoint_step_v307(state, bridge, hamiltonian, end_time,
                                               substeps=substeps, family_id=family_id, trajectory_id=trajectory_id, step=index,
                                               **metadata, basis_settings=basis_settings, midpoint_settings=midpoint_settings)
        trajectory.append(result)
        state, metadata = result.after, result.receipt["after_metadata"]
        if progress is not None:
            progress(result.receipt)
    return tuple(trajectory)


def serialize_checkpoint_trajectory_v307(trajectory):
    records = []
    for result in trajectory:
        record = dict(result.receipt)
        digest = record.pop("fingerprint")
        if record.get("schema") != SCHEMA_V307_CHECKPOINT or stable_hash_v294(record) != digest:
            raise ValueError("checkpoint receipt schema or fingerprint mismatch")
        if record["before"] != _state_payload(result.before) or record["after"] != _state_payload(result.after):
            raise ValueError("checkpoint states no longer match their receipt")
        if records:
            previous = records[-1]
            if previous["after"] != record["before"] or previous["after_metadata"] != record["before_metadata"]:
                raise ValueError("checkpoint states or metadata are not contiguous")
            if record["identity"] != {**previous["identity"], "step": previous["identity"]["step"]+1}:
                raise ValueError("checkpoint identities are not contiguous")
        records.append({**record, "fingerprint": digest})
    if not records:
        raise ValueError("checkpoint trajectory must not be empty")
    payload = {"schema": SCHEMA_V307_CHECKPOINT, "checkpoint_times": [r["checkpoint_time"] for r in records], "steps": records}
    return json.dumps({**payload, "fingerprint": stable_hash_v294(payload)}, sort_keys=True, separators=(",", ":"), allow_nan=False)+"\n"


def replay_checkpoint_trajectory_v307(serialized, bridge, hamiltonian):
    payload = json.loads(serialized)
    digest = payload.pop("fingerprint")
    if payload.get("schema") != SCHEMA_V307_CHECKPOINT or stable_hash_v294(payload) != digest:
        raise ValueError("checkpoint trajectory schema or fingerprint mismatch")
    records = payload["steps"]
    if not records or payload["checkpoint_times"] != [r["checkpoint_time"] for r in records]:
        raise ValueError("incomplete checkpoint time schedule")
    state, metadata = _state_from_payload(records[0]["before"]), records[0]["before_metadata"]
    first_identity = records[0]["identity"]
    replayed = []
    for offset, record in enumerate(records):
        rebuilt = molecular_checkpoint_step_v307(
            state, bridge, hamiltonian, record["checkpoint_time"], substeps=record["substeps"],
            **{**first_identity, "step": first_identity["step"]+offset}, **metadata,
            basis_settings=_restore_settings(record["basis_settings"]),
            midpoint_settings=MolecularMidpointSettingsV302(**record["midpoint_settings"]),
        )
        if rebuilt.receipt != record:
            raise ValueError("checkpoint replay does not reproduce stored receipt")
        replayed.append(rebuilt)
        state, metadata = rebuilt.after, rebuilt.receipt["after_metadata"]
    return tuple(replayed)
