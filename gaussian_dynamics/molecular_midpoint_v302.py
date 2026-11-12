"""Fixed-basis, full-parameter implicit midpoint and deterministic replay.

This propagates the v0.30.1 finite-grid projected model. It does not activate
spawning, change patches, renormalize states, or certify an ab-initio provider.
"""
from dataclasses import asdict, dataclass
import json

import numpy as np

from .correlated_gaussian_tdvp_v270 import (
    CorrelatedGaussianSpinorStateV270,
    pack_correlated_parameters_v270, state_from_correlated_parameters_v270,
)
from .molecular_frame_bridge_v300 import _hash_v300
from .molecular_hamiltonian_v301 import build_molecular_tdvp_snapshot_v301

SCHEMA_V302 = "gnd-molecular-midpoint-v0.30.2-dev"


@dataclass(frozen=True)
class MolecularMidpointSettingsV302:
    absolute_tolerance: float = 2e-10
    relative_tolerance: float = 2e-9
    tangent_step: float = 2e-6
    maximum_iterations: int = 24

    def validate(self):
        for name in ("absolute_tolerance", "relative_tolerance", "tangent_step"):
            value = getattr(self, name)
            if isinstance(value, (bool, np.bool_)) or not np.isfinite(value) or value <= 0:
                raise ValueError(f"{name} must be finite and positive")
        if type(self.maximum_iterations) is not int or self.maximum_iterations < 1:
            raise ValueError("maximum_iterations must be a positive integer")
        return self


def _restore(theta, template, time):
    return state_from_correlated_parameters_v270(
        theta, ngaussian=template.ngaussian, ndim=template.ndim,
        nstate=template.nstate, time_au=time,
    )


def _state_payload(state):
    # Store the actual matrices: log/exp of a general SPD width is not a
    # bit-exact round trip, even in the same runtime.
    return {"q": state.q.tolist(), "p": state.p.tolist(),
            "width_matrices": state.width_matrices.tolist(),
            "chirp_matrices": state.chirp_matrices.tolist(),
            "coefficients_real": state.coefficients.real.tolist(),
            "coefficients_imag": state.coefficients.imag.tolist(),
            "time_au": float(state.time_au)}


def _state_from_payload(payload):
    return CorrelatedGaussianSpinorStateV270(
        payload["q"], payload["p"], payload["width_matrices"], payload["chirp_matrices"],
        np.array(payload["coefficients_real"]) + 1j*np.array(payload["coefficients_imag"]),
        time_au=payload["time_au"],
    ).validate(require_normalized=False)


def _observables(evaluation, quadrature):
    psi, hpsi = evaluation.psi, evaluation.hpsi
    norm = float(np.sum(quadrature.weights[:, None] * np.abs(psi)**2))
    if not np.isfinite(norm) or norm <= 0:
        raise ValueError("finite-grid norm must be positive")
    energy = float(np.sum(quadrature.weights[:, None] * psi.conj() * hpsi).real / norm)
    return {"norm": norm, "energy_hartree": energy}


@dataclass(frozen=True)
class MolecularMidpointStepV302:
    before: object
    after: object
    receipt: dict


def molecular_midpoint_step_v302(state, bridge, hamiltonian, dt, *,
                                 settings=MolecularMidpointSettingsV302()):
    settings.validate()
    if isinstance(dt, (bool, np.bool_)) or not np.isfinite(dt) or dt == 0:
        raise ValueError("dt must be finite and nonzero")
    state.validate(require_normalized=False)
    start = pack_correlated_parameters_v270(state)
    evaluate = lambda s: build_molecular_tdvp_snapshot_v301(
        s, bridge, hamiltonian, relative_step=settings.tangent_step,
    )
    initial = evaluate(state)
    trial = start + dt * initial.snapshot.velocity
    end_time = float(state.time_au + dt)
    if not np.isfinite(end_time) or end_time == state.time_au:
        raise ValueError("time increment is not representable")
    for iteration in range(1, settings.maximum_iterations + 1):
        midpoint = _restore((start + trial)/2, state, state.time_au + dt/2)
        evaluation = evaluate(midpoint)
        residual = trial - start - dt * evaluation.snapshot.velocity
        scale = settings.absolute_tolerance + settings.relative_tolerance * np.maximum(abs(start), abs(trial))
        scaled_residual = float(np.max(abs(residual)/scale))
        if np.isfinite(scaled_residual) and scaled_residual <= 1.0:
            break
        trial = start + dt * evaluation.snapshot.velocity
    else:
        raise ValueError("implicit midpoint failed to converge; reduce dt")

    final_state = _restore(trial, state, end_time)
    final = evaluate(final_state)
    before_obs = _observables(initial, hamiltonian.quadrature)
    after_obs = _observables(final, hamiltonian.quadrature)
    audits = []
    for q0, q1 in zip(state.q, final_state.q):
        if np.array_equal(q0, q1):
            audits.append({"kind": "identity", "start": q0.tolist(), "end": q1.tolist()})
        else:
            audits.append({"kind": "audited-segment", "audit": bridge.audit_segment(q0, q1).as_dict()})
    solve = evaluation.snapshot.solve_receipt
    receipt = {
        "schema": SCHEMA_V302, "dt": float(dt), "iterations": iteration,
        "hamiltonian_fingerprint": hamiltonian.fingerprint,
        "bridge_fingerprint": bridge.bridge_fingerprint,
        "settings": asdict(settings),
        "before": _state_payload(state), "after": _state_payload(final_state),
        "midpoint_residual_inf": float(np.max(abs(residual))),
        "scaled_midpoint_residual": scaled_residual,
        "midpoint_metric_rank": solve.rank,
        "midpoint_metric_condition": solve.condition_number,
        "midpoint_svd_residual": solve.residual_relative,
        "midpoint_null_rhs_relative": solve.null_rhs_relative,
        "midpoint_tdse_defect_norm": float(np.linalg.norm(evaluation.snapshot.residual_weighted)),
        "before_observables": before_obs, "after_observables": after_obs,
        "center_transport_audits": audits,
    }
    receipt["fingerprint"] = _hash_v300(receipt)
    return MolecularMidpointStepV302(state, final_state, receipt)


def run_molecular_midpoint_v302(state, bridge, hamiltonian, *, dt, steps,
                               settings=MolecularMidpointSettingsV302()):
    if type(steps) is not int or steps < 1:
        raise ValueError("steps must be a positive integer")
    trajectory = []
    for _ in range(steps):
        step = molecular_midpoint_step_v302(state, bridge, hamiltonian, dt, settings=settings)
        trajectory.append(step)
        state = step.after
    return tuple(trajectory)


def serialize_molecular_trajectory_v302(trajectory):
    receipts = []
    for step in trajectory:
        receipt = dict(step.receipt)
        fingerprint = receipt.pop("fingerprint")
        if _hash_v300(receipt) != fingerprint:
            raise ValueError("trajectory receipt fingerprint mismatch")
        if _state_payload(step.before) != receipt["before"] or _state_payload(step.after) != receipt["after"]:
            raise ValueError("trajectory state no longer matches receipt")
        if receipts and receipts[-1]["after"] != receipt["before"]:
            raise ValueError("trajectory steps are not contiguous")
        receipt["fingerprint"] = fingerprint
        receipts.append(receipt)
    if not receipts:
        raise ValueError("trajectory must not be empty")
    payload = {"schema": SCHEMA_V302, "steps": receipts}
    return json.dumps({**payload, "fingerprint": _hash_v300(payload)}, sort_keys=True,
                      separators=(",", ":"), allow_nan=False) + "\n"


def replay_molecular_trajectory_v302(serialized, bridge, hamiltonian):
    """Rerun every solve and audit; exact receipts require the same numeric runtime.

    Fingerprints detect accidental changes, not malicious re-signing. Replayed
    dynamics, residuals and provider fingerprints are independently recomputed.
    """
    payload = json.loads(serialized)
    digest = payload.pop("fingerprint")
    if payload.get("schema") != SCHEMA_V302 or _hash_v300(payload) != digest:
        raise ValueError("trajectory schema or fingerprint mismatch")
    records = payload["steps"]
    if not records:
        raise ValueError("trajectory must not be empty")
    state = _state_from_payload(records[0]["before"])
    replayed = []
    for record in records:
        if record["hamiltonian_fingerprint"] != hamiltonian.fingerprint:
            raise ValueError("replay Hamiltonian fingerprint mismatch")
        step = molecular_midpoint_step_v302(
            state, bridge, hamiltonian, record["dt"],
            settings=MolecularMidpointSettingsV302(**record["settings"]),
        )
        if step.receipt != record:
            raise ValueError("trajectory replay does not reproduce stored receipt")
        replayed.append(step)
        state = step.after
    return tuple(replayed)
