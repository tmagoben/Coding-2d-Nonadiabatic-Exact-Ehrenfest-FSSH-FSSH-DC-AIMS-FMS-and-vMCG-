"""Implicit midpoint with dormant shapes and audited newborn activation.

All electronic coefficients remain active. Dormant q, p, width and chirp arrays
are copied exactly through each parameter restoration, avoiding log/exp drift.
"""
from dataclasses import asdict, dataclass, replace
import json

import numpy as np

from .correlated_gaussian_tdvp_v270 import pack_correlated_parameters_v270
from .curved_basis_adaptation_v292 import CurvedBasisSettingsV292
from .curved_basis_adaptation_v293 import active_parameter_indices_v293
from .molecular_hamiltonian_v301 import build_molecular_tdvp_snapshot_v301
from .molecular_midpoint_v302 import (
    MolecularMidpointSettingsV302, _restore, _state_payload, _state_from_payload, _observables,
)
from .pre_ml_contract_v294 import stable_hash_v294

SCHEMA_V306_ACTIVE = "gnd-active-molecular-midpoint-v0.30.6-dev"
SCHEMA_V306_ACTIVATION = "gnd-molecular-shape-activation-v0.30.6-dev"


def _mask(state, active_shape_mask):
    mask = np.ones(state.ngaussian, dtype=bool) if active_shape_mask is None else np.asarray(active_shape_mask)
    active_parameter_indices_v293(state, mask)
    return mask.copy()


def _restore_active(theta, state, time, mask):
    restored = _restore(theta, state, time)
    values = {}
    for name in ("q", "p", "width_matrices", "chirp_matrices"):
        value = getattr(restored, name).copy()
        value[~mask] = getattr(state, name)[~mask]
        values[name] = value
    return replace(restored, **values)


@dataclass(frozen=True)
class MolecularActiveStepV306:
    before: object
    after: object
    receipt: dict


def molecular_active_midpoint_step_v306(state, bridge, hamiltonian, dt, *,
                                       active_shape_mask=None, settings=MolecularMidpointSettingsV302()):
    """Solve only active parameter equations; no renormalization or topology change."""
    state.validate(require_normalized=False)
    settings.validate()
    if isinstance(dt, (bool, np.bool_)) or not np.isfinite(dt) or dt == 0:
        raise ValueError("dt must be finite and nonzero")
    mask = _mask(state, active_shape_mask)
    indices = active_parameter_indices_v293(state, mask)
    start = pack_correlated_parameters_v270(state)
    def evaluate(current):
        return build_molecular_tdvp_snapshot_v301(
            current, bridge, hamiltonian, relative_step=settings.tangent_step, active_shape_mask=mask,
        )
    initial = evaluate(state)
    trial = start.copy()
    trial[indices] += dt * initial.snapshot.velocity
    end_time = float(state.time_au + dt)
    if not np.isfinite(end_time) or end_time == state.time_au:
        raise ValueError("time increment is not representable")
    for iteration in range(1, settings.maximum_iterations + 1):
        midpoint = _restore_active((start + trial)/2, state, state.time_au + dt/2, mask)
        evaluation = evaluate(midpoint)
        residual = trial[indices] - start[indices] - dt * evaluation.snapshot.velocity
        scale = settings.absolute_tolerance + settings.relative_tolerance * np.maximum(abs(start[indices]), abs(trial[indices]))
        scaled = float(np.max(abs(residual)/scale))
        if np.isfinite(scaled) and scaled <= 1.0:
            break
        trial[indices] = start[indices] + dt * evaluation.snapshot.velocity
    else:
        raise ValueError("active implicit midpoint failed to converge; reduce dt")
    after = _restore_active(trial, state, end_time, mask)
    final = evaluate(after)
    audits = []
    for q0, q1 in zip(state.q, after.q):
        if np.array_equal(q0, q1):
            audits.append({"kind": "identity", "start": q0.tolist(), "end": q1.tolist()})
        else:
            audits.append({"kind": "audited-segment", "audit": bridge.audit_segment(q0, q1).as_dict()})
    solve = evaluation.snapshot.solve_receipt
    receipt = {
        "schema": SCHEMA_V306_ACTIVE, "dt": float(dt), "iterations": iteration,
        "active_shape_mask": mask.tolist(), "active_parameter_indices": indices.tolist(),
        "hamiltonian_fingerprint": hamiltonian.fingerprint, "bridge_fingerprint": bridge.bridge_fingerprint,
        "settings": asdict(settings), "before": _state_payload(state), "after": _state_payload(after),
        "midpoint_residual_inf": float(np.max(abs(residual))), "scaled_midpoint_residual": scaled,
        "midpoint_metric_rank": solve.rank, "midpoint_metric_condition": solve.condition_number,
        "midpoint_svd_residual": solve.residual_relative, "midpoint_null_rhs_relative": solve.null_rhs_relative,
        "midpoint_tdse_defect_norm": float(np.linalg.norm(evaluation.snapshot.residual_weighted)),
        "before_observables": _observables(initial, hamiltonian.quadrature),
        "after_observables": _observables(final, hamiltonian.quadrature), "center_transport_audits": audits,
    }
    receipt["fingerprint"] = stable_hash_v294(receipt)
    return MolecularActiveStepV306(state, after, receipt)


def serialize_active_step_v306(step):
    receipt = dict(step.receipt)
    digest = receipt.pop("fingerprint")
    if receipt.get("schema") != SCHEMA_V306_ACTIVE or stable_hash_v294(receipt) != digest:
        raise ValueError("active step schema or fingerprint mismatch")
    if _state_payload(step.before) != receipt["before"] or _state_payload(step.after) != receipt["after"]:
        raise ValueError("active step state no longer matches receipt")
    return json.dumps({**receipt, "fingerprint": digest}, sort_keys=True, separators=(",", ":"), allow_nan=False) + "\n"


def replay_active_step_v306(serialized, bridge, hamiltonian):
    receipt = json.loads(serialized)
    payload = dict(receipt)
    digest = payload.pop("fingerprint")
    if payload.get("schema") != SCHEMA_V306_ACTIVE or stable_hash_v294(payload) != digest:
        raise ValueError("active step schema or fingerprint mismatch")
    if receipt["hamiltonian_fingerprint"] != hamiltonian.fingerprint:
        raise ValueError("replay Hamiltonian fingerprint mismatch")
    result = molecular_active_midpoint_step_v306(
        _state_from_payload(receipt["before"]), bridge, hamiltonian, receipt["dt"],
        active_shape_mask=receipt["active_shape_mask"], settings=MolecularMidpointSettingsV302(**receipt["settings"]),
    )
    if result.receipt != receipt:
        raise ValueError("active midpoint replay does not reproduce receipt")
    return result


def activate_molecular_shapes_v306(state, bridge, hamiltonian, active_shape_mask, *,
                                  packet_ids=None, settings=CurvedBasisSettingsV292(), tangent_step=2e-6):
    """Audit dormant blocks in stable-ID order without modifying the state.

    Condition is measured on new shape tangents after projecting out the current
    tangent. Velocity amplification uses shape-coordinate speeds, excluding the
    representation-dependent electronic coefficient derivatives.
    """
    state.validate(require_normalized=False)
    settings.validate()
    mask = _mask(state, active_shape_mask)
    if isinstance(tangent_step, (bool, np.bool_)) or not np.isfinite(tangent_step) or tangent_step <= 0:
        raise ValueError("tangent_step must be finite and positive")
    if bridge.bridge_fingerprint != hamiltonian.bridge_fingerprint:
        raise ValueError("Hamiltonian and molecular bridge fingerprints disagree")
    before = mask.copy()
    ids = tuple(f"g{i:06d}" for i in range(state.ngaussian)) if packet_ids is None else tuple(packet_ids)
    if len(ids) != state.ngaussian or any(not isinstance(x, str) or not x.strip() or x != x.strip() for x in ids) or len(set(ids)) != len(ids):
        raise ValueError("packet IDs must be nonempty, unique and complete")
    def evaluate(current_mask):
        return build_molecular_tdvp_snapshot_v301(state, bridge, hamiltonian,
                                                relative_step=tangent_step, active_shape_mask=current_mask).snapshot
    base = evaluate(mask) if not np.all(mask) else None
    trials = []
    for i in sorted(np.flatnonzero(~mask), key=lambda k: ids[k]):
        population = float(np.vdot(state.coefficients[i], state.coefficients[i]).real)
        row = {"packet_id": ids[i], "packet_index": int(i), "coefficient_norm_squared": population,
               "activated": False, "reason": None, "projected_shape_rank": None,
               "projected_shape_condition": None, "shape_velocity_amplification": None,
               "projected_shape_metric_spectrum": None, "projected_shape_cutoff": None,
               "base_shape_speed": None, "trial_shape_speed": None,
               "trial_metric_rank": None, "trial_svd_residual": None}
        if population < settings.inherited.newborn_activation_population:
            row["reason"] = "insufficient_coefficient_population"
            trials.append(row)
            continue
        trial_mask = mask.copy()
        trial_mask[i] = True
        try:
            trial = evaluate(trial_mask)
        except ValueError as error:
            row.update(reason="incompatible_trial_metric", detail=str(error))
            trials.append(row)
            continue
        base_indices = active_parameter_indices_v293(state, mask)
        trial_indices = active_parameter_indices_v293(state, trial_mask)
        novel_columns = np.flatnonzero(~np.isin(trial_indices, base_indices))
        j = np.vstack([base.jacobian.real, base.jacobian.imag])
        k = np.vstack([trial.jacobian[:, novel_columns].real, trial.jacobian[:, novel_columns].imag])
        u, s, _ = np.linalg.svd(j, full_matrices=False)
        keep = s*s > max(settings.tangent_absolute_cutoff, settings.tangent_relative_cutoff * float(s[0]**2))
        novel = k - u[:, keep] @ (u[:, keep].T @ k)
        spectrum = np.linalg.svd(novel, compute_uv=False)**2
        cutoff = max(settings.tangent_absolute_cutoff, settings.tangent_relative_cutoff * float(spectrum[0]))
        rank = int(np.count_nonzero(spectrum > cutoff))
        condition = float(spectrum[0]/spectrum[rank-1]) if rank else None
        first_shape = 2 * state.ngaussian * state.nstate
        base_speed = float(np.linalg.norm(base.velocity[base_indices >= first_shape]))
        trial_speed = float(np.linalg.norm(trial.velocity[trial_indices >= first_shape]))
        amplification = trial_speed / max(base_speed, 1e-30)
        row.update(projected_shape_rank=rank, projected_shape_condition=condition,
                   projected_shape_metric_spectrum=spectrum.tolist(), projected_shape_cutoff=cutoff,
                   base_shape_speed=base_speed, trial_shape_speed=trial_speed,
                   shape_velocity_amplification=amplification if np.isfinite(amplification) else None,
                   trial_metric_rank=trial.solve_receipt.rank, trial_svd_residual=trial.solve_receipt.residual_relative)
        if rank != len(novel_columns):
            row["reason"] = "incomplete_new_shape_tangent"
        elif condition > settings.inherited.maximum_activation_condition_number:
            row["reason"] = "ill_conditioned_new_shape_tangent"
        elif not np.isfinite(amplification) or amplification > settings.inherited.maximum_activation_velocity_amplification:
            row["reason"] = "shape_velocity_amplification"
        else:
            row.update(activated=True, reason="activation_gates_passed")
            mask, base = trial_mask, trial
        trials.append(row)
    receipt = {"schema": SCHEMA_V306_ACTIVATION, "state_fingerprint": stable_hash_v294(_state_payload(state)),
               "hamiltonian_fingerprint": hamiltonian.fingerprint, "bridge_fingerprint": bridge.bridge_fingerprint,
               "settings": asdict(settings), "tangent_step": float(tangent_step), "packet_ids": list(ids),
               "mask_before": before.tolist(), "mask_after": mask.tolist(), "trials": trials}
    receipt["fingerprint"] = stable_hash_v294(receipt)
    return mask, receipt
