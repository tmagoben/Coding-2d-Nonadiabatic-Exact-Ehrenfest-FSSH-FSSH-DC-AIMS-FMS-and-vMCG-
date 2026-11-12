"""Audited propagation -> ageing -> activation -> scheduled topology checkpoints.

Topology is deferred while any shape block remains dormant. This explicitly
versioned policy avoids scoring a fully active tangent for a constrained state.
"""
from dataclasses import asdict, dataclass
import json

import numpy as np

from .curved_basis_adaptation_v292 import CurvedBasisSettingsV292
from .molecular_active_v306 import molecular_active_midpoint_step_v306, activate_molecular_shapes_v306
from .molecular_admission_v305 import _metadata, evaluate_molecular_admission_v305
from .molecular_midpoint_v302 import MolecularMidpointSettingsV302, _state_payload, _state_from_payload
from .molecular_spawn_v304 import _name, _restore_settings
from .pre_ml_contract_v294 import stable_hash_v294

SCHEMA_V306_ADAPTIVE = "gnd-molecular-adaptive-trajectory-v0.30.6-dev"
POLICY_V306 = "propagate; age; activate; scheduled merge/prune/spawn only when all shapes active"


@dataclass(frozen=True)
class MolecularAdaptiveStepV306:
    before: object
    after: object
    receipt: dict


def molecular_adaptive_step_v306(
    state, bridge, hamiltonian, dt, *, family_id, trajectory_id, step,
    packet_ids=None, packet_ages=None, active_shape_mask=None, next_packet_serial=None,
    basis_settings=CurvedBasisSettingsV292(), midpoint_settings=MolecularMidpointSettingsV302(),
):
    """One forward macro-step. Failure returns no accepted event or mutated state."""
    state.validate(require_normalized=False)
    basis_settings.validate()
    midpoint_settings.validate()
    if isinstance(dt, (bool, np.bool_)) or not np.isfinite(dt) or dt <= 0:
        raise ValueError("adaptive dt must be finite and positive")
    if type(step) is not int or step < 1:
        raise ValueError("adaptive step must be a positive integer")
    family_id, trajectory_id = _name(family_id, "family_id"), _name(trajectory_id, "trajectory_id")
    ids, ages, mask, serial = _metadata(state, packet_ids, packet_ages, active_shape_mask, next_packet_serial)
    before_metadata = dict(packet_ids=list(ids), packet_ages=list(ages), active_shape_mask=list(mask), next_packet_serial=serial)
    propagation = molecular_active_midpoint_step_v306(
        state, bridge, hamiltonian, dt, active_shape_mask=mask, settings=midpoint_settings,
    )
    active, activation = activate_molecular_shapes_v306(
        propagation.after, bridge, hamiltonian, mask, packet_ids=ids,
        settings=basis_settings, tangent_step=midpoint_settings.tangent_step,
    )
    after_metadata = dict(packet_ids=list(ids), packet_ages=[age + 1 for age in ages],
                          active_shape_mask=active.tolist(), next_packet_serial=serial)
    after = propagation.after
    topology = None
    scheduled = step % basis_settings.inherited.adapt_every_steps == 0
    if not scheduled:
        status = "not_scheduled"
    elif not np.all(active):
        status = "deferred_dormant_shapes"
    else:
        status = "evaluated"
        admission = evaluate_molecular_admission_v305(
            after, bridge, hamiltonian, family_id=family_id, trajectory_id=trajectory_id,
            step=step, **after_metadata, settings=basis_settings, tangent_step=midpoint_settings.tangent_step,
        )
        after, topology = admission.after, admission.receipt
        after_metadata = topology["after_metadata"]
    receipt = {
        "schema": SCHEMA_V306_ADAPTIVE, "policy": POLICY_V306,
        "identity": dict(family_id=family_id, trajectory_id=trajectory_id, step=step),
        "dt": float(dt), "basis_settings": asdict(basis_settings), "midpoint_settings": asdict(midpoint_settings),
        "before": _state_payload(state), "before_metadata": before_metadata,
        "propagation": propagation.receipt, "activation": activation,
        "topology_status": status, "topology": topology,
        "after": _state_payload(after), "after_metadata": after_metadata,
        "dataset_certified": False, "training_eligible": False,
    }
    receipt["fingerprint"] = stable_hash_v294(receipt)
    return MolecularAdaptiveStepV306(state, after, receipt)


def run_molecular_adaptive_v306(state, bridge, hamiltonian, *, time_steps,
                               family_id, trajectory_id, initial_step=0,
                               packet_ids=None, packet_ages=None, active_shape_mask=None, next_packet_serial=None,
                               basis_settings=CurvedBasisSettingsV292(), midpoint_settings=MolecularMidpointSettingsV302(),
                               progress=None):
    """Advance an explicit forward schedule; no hidden retries or time-step changes."""
    times = tuple(time_steps)
    if not times or any(isinstance(dt, (bool, np.bool_)) or not np.isfinite(dt) or dt <= 0 for dt in times):
        raise ValueError("time_steps must contain finite positive increments")
    if type(initial_step) is not int or initial_step < 0:
        raise ValueError("initial_step must be a nonnegative integer")
    metadata = dict(packet_ids=packet_ids, packet_ages=packet_ages,
                    active_shape_mask=active_shape_mask, next_packet_serial=next_packet_serial)
    trajectory = []
    for offset, dt in enumerate(times, 1):
        accepted = molecular_adaptive_step_v306(
            state, bridge, hamiltonian, dt, family_id=family_id, trajectory_id=trajectory_id,
            step=initial_step + offset, **metadata, basis_settings=basis_settings, midpoint_settings=midpoint_settings,
        )
        trajectory.append(accepted)
        state, metadata = accepted.after, accepted.receipt["after_metadata"]
        if progress is not None:
            event = accepted.receipt["topology"]
            progress({"step": initial_step + offset, "time_au": state.time_au, "packet_count": state.ngaussian,
                      "active_shape_mask": metadata["active_shape_mask"],
                      "event": event["event_kind"] if event is not None else accepted.receipt["topology_status"],
                      "activation_trials": accepted.receipt["activation"]["trials"]})
    return tuple(trajectory)


def serialize_adaptive_trajectory_v306(trajectory):
    records = []
    for accepted in trajectory:
        receipt = dict(accepted.receipt)
        digest = receipt.pop("fingerprint")
        if receipt.get("schema") != SCHEMA_V306_ADAPTIVE or stable_hash_v294(receipt) != digest:
            raise ValueError("adaptive receipt schema or fingerprint mismatch")
        if _state_payload(accepted.before) != receipt["before"] or _state_payload(accepted.after) != receipt["after"]:
            raise ValueError("adaptive trajectory state no longer matches receipt")
        if records:
            previous = records[-1]
            if previous["after"] != receipt["before"] or previous["after_metadata"] != receipt["before_metadata"]:
                raise ValueError("adaptive states or metadata are not contiguous")
            expected = {**previous["identity"], "step": previous["identity"]["step"] + 1}
            if receipt["identity"] != expected:
                raise ValueError("adaptive step identities are not contiguous")
        records.append({**receipt, "fingerprint": digest})
    if not records:
        raise ValueError("adaptive trajectory must not be empty")
    payload = {"schema": SCHEMA_V306_ADAPTIVE, "time_steps": [r["dt"] for r in records], "steps": records}
    return json.dumps({**payload, "fingerprint": stable_hash_v294(payload)}, sort_keys=True, separators=(",", ":"), allow_nan=False) + "\n"


def replay_adaptive_trajectory_v306(serialized, bridge, hamiltonian):
    payload = json.loads(serialized)
    digest = payload.pop("fingerprint")
    if payload.get("schema") != SCHEMA_V306_ADAPTIVE or stable_hash_v294(payload) != digest:
        raise ValueError("adaptive trajectory schema or fingerprint mismatch")
    records = payload["steps"]
    if not records or payload["time_steps"] != [r["dt"] for r in records]:
        raise ValueError("adaptive trajectory has an incomplete time schedule")
    state = _state_from_payload(records[0]["before"])
    metadata = records[0]["before_metadata"]
    first_identity = records[0]["identity"]
    replayed = []
    for offset, record in enumerate(records):
        identity = {**first_identity, "step": first_identity["step"] + offset}
        rebuilt = molecular_adaptive_step_v306(
            state, bridge, hamiltonian, record["dt"], **identity, **metadata,
            basis_settings=_restore_settings(record["basis_settings"]),
            midpoint_settings=MolecularMidpointSettingsV302(**record["midpoint_settings"]),
        )
        if rebuilt.receipt != record:
            raise ValueError("adaptive replay does not reproduce stored step")
        replayed.append(rebuilt)
        state, metadata = rebuilt.after, rebuilt.receipt["after_metadata"]
    return tuple(replayed)
