"""Audited, single-event molecular topology decisions on the fixed grid.

Residual screening and final admission are separate outcomes. A newborn has
active coefficients and dormant shape parameters; this module does not advance
time, activate newborn shapes, or certify an adaptive trajectory.
"""
from dataclasses import asdict, dataclass
import json

import numpy as np

from .correlated_basis_adaptation_v270 import MultidimensionalSpawnCandidateV270
from .correlated_gaussian_tdvp_v270 import CorrelatedGaussianSpinorStateV270
from .curved_basis_adaptation_v292 import CurvedBasisSettingsV292, _canon
from .curved_basis_adaptation_v293 import (
    evaluate_curved_state_v293, packet_subspace_overlap_v293,
    project_curved_state_v293,
)
from .molecular_midpoint_v302 import _state_payload, _state_from_payload
from .molecular_spawn_v304 import (
    _checkpoint_connection, _restore_settings,
    evaluate_molecular_spawn_checkpoint_v304,
)
from .pre_ml_contract_v294 import stable_hash_v294

SCHEMA_V305 = "gnd-molecular-admission-v0.30.5-dev"
POLICY_V305 = "merge-prune-spawn; one event; selected candidate only; projection and energy gates"


@dataclass(frozen=True)
class MolecularAdmissionV305:
    after: CorrelatedGaussianSpinorStateV270
    receipt: dict


def _metadata(state, packet_ids, packet_ages, active_shape_mask, next_packet_serial):
    n = state.ngaussian
    ids = tuple(f"g{i:06d}" for i in range(n)) if packet_ids is None else tuple(packet_ids)
    ages = (0,) * n if packet_ages is None else tuple(packet_ages)
    mask = (True,) * n if active_shape_mask is None else tuple(active_shape_mask)
    serial = n if next_packet_serial is None else next_packet_serial
    if len(ids) != n or any(not isinstance(x, str) or not x.strip() or x != x.strip() for x in ids) or len(set(ids)) != n:
        raise ValueError("packet IDs must be nonempty, unique and complete")
    if len(ages) != n or any(type(x) is not int or x < 0 for x in ages):
        raise ValueError("packet ages must be nonnegative integers")
    if len(mask) != n or any(not isinstance(x, (bool, np.bool_)) for x in mask):
        raise ValueError("active shape mask must contain one Boolean per packet")
    if type(serial) is not int or serial < 0 or f"g{serial:06d}" in ids:
        raise ValueError("next packet serial must be nonnegative and unused")
    return ids, ages, tuple(bool(x) for x in mask), serial


def _energy(state, connection, hamiltonian):
    psi = evaluate_curved_state_v293(state, connection, hamiltonian.quadrature)
    weighted = hamiltonian.quadrature.weights[:, None]
    norm = float(np.sum(weighted * abs(psi)**2))
    if not np.isfinite(norm) or norm <= 0:
        raise ValueError("nonpositive or nonfinite grid norm")
    value = np.sum(weighted * psi.conj() * hamiltonian.apply(psi)) / norm
    if not np.isfinite(value) or abs(value.imag) > 2e-10 * max(1.0, abs(value.real)):
        raise ValueError("grid energy must be finite and real")
    return float(value.real)


def _project_event(state, target, connection, hamiltonian, settings, source_energy):
    """Return explicit rejection evidence; an energy evaluation error is fatal."""
    try:
        after, projection = project_curved_state_v293(
            state, target, connection, hamiltonian.quadrature, settings=settings,
        )
    except ValueError as error:
        return None, {"outcome": "projection_rejected", "detail": str(error),
                      "projection": None, "energy_before_hartree": source_energy,
                      "energy_after_hartree": None, "energy_jump_hartree": None}
    target_energy = _energy(after, connection, hamiltonian)
    delta = target_energy - source_energy
    loss = projection.relative_loss
    # A nonfinite receipt must never satisfy a gate through NaN comparisons.
    if not all(np.isfinite(x) for x in (loss, projection.fidelity, projection.condition_number)):
        raise ValueError("nonfinite projection receipt")
    outcome = ("projection_rejected" if loss > settings.inherited.maximum_projection_loss
               else "energy_rejected" if abs(delta) > settings.inherited.maximum_event_energy_jump_hartree
               else "admitted")
    evidence = {"outcome": outcome, "detail": None, "projection": _canon(asdict(projection)),
                "energy_before_hartree": source_energy, "energy_after_hartree": target_energy,
                "energy_jump_hartree": delta}
    return after if outcome == "admitted" else None, evidence


def evaluate_molecular_admission_v305(
    state, bridge, hamiltonian, *, family_id, trajectory_id, step,
    packet_ids=None, packet_ages=None, active_shape_mask=None, next_packet_serial=None,
    settings=CurvedBasisSettingsV292(), tangent_step=2e-6,
):
    """Evaluate all residual candidates, then apply one ordered topology event.

    This is an explicitly requested adaptation checkpoint, not a timestep
    scheduler. Ages are supplied by the caller and are not incremented here.
    Nonselected candidates have not failed projection or energy gates.
    """
    state.validate(require_normalized=False)
    settings.validate()
    ids, ages, mask, serial = _metadata(state, packet_ids, packet_ages, active_shape_mask, next_packet_serial)
    checkpoint = evaluate_molecular_spawn_checkpoint_v304(
        state, bridge, hamiltonian, family_id=family_id, trajectory_id=trajectory_id,
        step=step, packet_ids=ids, settings=settings, tangent_step=tangent_step,
    )
    connection = _checkpoint_connection(bridge)
    source_energy = _energy(state, connection, hamiltonian)
    inherited = settings.inherited
    rows = [{"serial": row["serial"], "residual_accepted": row["record"]["accepted_by_residual_oracle"],
             "final_admitted": False, "outcome": "not_selected" if row["record"]["accepted_by_residual_oracle"] else "residual_rejected",
             "attempt": None} for row in checkpoint["candidates"]]
    attempts = []
    after = state
    after_ids, after_ages, after_mask = ids, ages, mask
    event_kind, event_reason = "none", "no_residual_candidate"
    selected_serial = added_id = removed_id = None

    def removal(remove, kind):
        keep = [i for i in range(state.ngaussian) if i != remove]
        target = CorrelatedGaussianSpinorStateV270(
            state.q[keep], state.p[keep], state.width_matrices[keep],
            state.chirp_matrices[keep], state.coefficients[keep], state.time_au,
        )
        projected, audit = _project_event(state, target, connection, hamiltonian, settings, source_energy)
        attempts.append({"kind": kind, "removed_packet_id": ids[remove], **audit})
        return keep, projected

    # Stable IDs make removal order independent of input packet order.
    ordered = sorted(range(state.ngaussian), key=lambda i: ids[i])
    if inherited.merge_enabled:
        for position, i in enumerate(ordered):
            for j in ordered[position + 1:]:
                if packet_subspace_overlap_v293(state, i, j, connection, hamiltonian.quadrature) < inherited.minimum_merge_overlap:
                    continue
                # On equal populations retain the lexicographically smaller ID.
                pi, pj = (float(np.vdot(state.coefficients[k], state.coefficients[k]).real) for k in (i, j))
                remove = j if pj <= pi else i
                keep, projected = removal(remove, "merge")
                if projected is not None:
                    after, event_kind, removed_id = projected, "merge", ids[remove]
                    after_ids, after_ages, after_mask = (tuple(values[k] for k in keep) for values in (ids, ages, mask))
                    break
            if event_kind != "none":
                break
    if event_kind == "none" and inherited.prune_enabled and state.ngaussian > 1:
        order = sorted(ordered, key=lambda i: (float(np.sum(abs(state.coefficients[i])**2)), ids[i]))
        for remove in order:
            if ages[remove] < inherited.minimum_prune_age or float(np.sum(abs(state.coefficients[remove])**2)) > inherited.prune_population_threshold:
                continue
            keep, projected = removal(remove, "prune")
            if projected is not None:
                after, event_kind, removed_id = projected, "prune", ids[remove]
                after_ids, after_ages, after_mask = (tuple(values[k] for k in keep) for values in (ids, ages, mask))
                break

    blocker = ("earlier_topology_event" if event_kind != "none" else
               "spawn_disabled" if not inherited.spawn_enabled else
               "dormant_shape_pending" if not all(mask) else
               "maximum_packet_count" if state.ngaussian >= inherited.maximum_packet_count else None)
    if blocker is not None:
        event_reason = blocker
        for row in rows:
            if row["residual_accepted"]:
                row["outcome"] = blocker
    else:
        eligible = [row for row in checkpoint["candidates"] if row["record"]["accepted_by_residual_oracle"]]
        if eligible:
            chosen = min(eligible, key=lambda row: (
                -row["record"]["predicted_residual_reduction"], row["source_packet_id"],
                MultidimensionalSpawnCandidateV270(**row["candidate"]).canonical_key()[1:],
            ))
            selected_serial = chosen["serial"]
            candidate = MultidimensionalSpawnCandidateV270(**chosen["candidate"])
            target = CorrelatedGaussianSpinorStateV270(
                np.vstack([state.q, candidate.q]), np.vstack([state.p, candidate.p]),
                np.concatenate([state.width_matrices, candidate.width_matrices[None]]),
                np.concatenate([state.chirp_matrices, candidate.chirp_matrices[None]]),
                np.vstack([state.coefficients, np.zeros((1, state.nstate), complex)]), state.time_au,
            )
            projected, audit = _project_event(state, target, connection, hamiltonian, settings, source_energy)
            rows[selected_serial].update(outcome=audit["outcome"], attempt=audit, final_admitted=projected is not None)
            event_reason = audit["outcome"]
            if projected is not None:
                after, event_kind, added_id = projected, "spawn", f"g{serial:06d}"
                after_ids, after_ages, after_mask = ids + (added_id,), ages + (0,), mask + (False,)

    receipt = {
        "schema": SCHEMA_V305, "policy": POLICY_V305,
        "screening": checkpoint, "candidate_space": "full_complex_electronic_spinor",
        "decision_stage": "final_topology_admission", "target_state": None,
        "before_metadata": {"packet_ids": list(ids), "packet_ages": list(ages), "active_shape_mask": list(mask), "next_packet_serial": serial},
        "after_metadata": {"packet_ids": list(after_ids), "packet_ages": list(after_ages), "active_shape_mask": list(after_mask), "next_packet_serial": serial + (event_kind == "spawn")},
        "event_kind": event_kind, "event_reason": event_reason,
        "selected_candidate_serial": selected_serial, "added_packet_id": added_id,
        "removed_packet_id": removed_id, "removal_attempts": attempts,
        "candidate_decisions": rows, "after_state": _state_payload(after),
        "dataset_certified": False, "adaptive_trajectory_validated": False,
    }
    receipt["fingerprint"] = stable_hash_v294(receipt)
    return MolecularAdmissionV305(after, receipt)


def serialize_molecular_admission_v305(receipt):
    payload = dict(receipt)
    digest = payload.pop("fingerprint")
    if payload.get("schema") != SCHEMA_V305 or stable_hash_v294(payload) != digest:
        raise ValueError("admission schema or fingerprint mismatch")
    return json.dumps(receipt, sort_keys=True, separators=(",", ":"), allow_nan=False) + "\n"


def replay_molecular_admission_v305(serialized, bridge, hamiltonian):
    """Recompute residual scores, all attempted events and the final state."""
    receipt = json.loads(serialized)
    serialize_molecular_admission_v305(receipt)
    screening = receipt["screening"]
    rebuilt = evaluate_molecular_admission_v305(
        _state_from_payload(screening["state"]), bridge, hamiltonian,
        **screening["identity"], **receipt["before_metadata"],
        settings=_restore_settings(screening["settings"]), tangent_step=screening["tangent_step"],
    )
    if rebuilt.receipt != receipt:
        raise ValueError("admission replay does not reproduce stored checkpoint")
    return rebuilt
