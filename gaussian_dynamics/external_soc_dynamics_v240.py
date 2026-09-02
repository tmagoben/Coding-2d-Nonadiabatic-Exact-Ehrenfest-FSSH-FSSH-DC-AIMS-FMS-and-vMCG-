"""Admission-bound frozen-snapshot electronic dynamics for v0.24.0."""

from dataclasses import dataclass

import numpy as np
from scipy.linalg import expm

from .external_soc_admission_v240 import ExternalSOCAdmissionAuditV240
from .openmolcas_rassi_snapshot_v240 import ParsedOpenMolcasBundleV240


@dataclass(frozen=True)
class FrozenSnapshotCheckpointV240:
    bundle_fingerprint: str
    step: int
    time_au: float
    coefficients: np.ndarray


@dataclass(frozen=True)
class FrozenSnapshotDynamicsV240:
    bundle_fingerprint: str
    evidence_class: str
    soc_enabled: bool
    times_au: np.ndarray
    coefficients: np.ndarray
    norms: np.ndarray
    checkpoint: FrozenSnapshotCheckpointV240


def _validate_state_v240(vector, nstate):
    vector = np.asarray(vector, dtype=complex)
    if vector.shape != (nstate,) or not np.all(np.isfinite(vector)):
        raise ValueError("initial electronic state has incompatible or non-finite data.")
    norm = float(np.linalg.norm(vector))
    if norm == 0.0:
        raise ValueError("initial electronic state cannot be zero.")
    return vector / norm


def preview_frozen_snapshot_dynamics_v240(
    bundle,
    initial_coefficients,
    *,
    time_step_au,
    steps,
    soc_enabled=True,
    checkpoint=None,
):
    """Propagate a parsed snapshot and retain its evidence classification.

    This preview is usable for protocol fixtures, but its result is explicitly marked
    as such.  Production entry is through ``run_admitted_external_soc_dynamics_v240``.
    """

    if type(bundle) is not ParsedOpenMolcasBundleV240:
        raise TypeError("frozen-snapshot dynamics requires a parsed bundle.")
    if type(soc_enabled) is not bool:
        raise TypeError("soc_enabled must be a native Boolean.")
    dt = float(time_step_au)
    if not np.isfinite(dt) or dt <= 0.0:
        raise ValueError("time_step_au must be finite and positive.")
    if int(steps) != steps or steps < 0:
        raise ValueError("steps must be a nonnegative integer.")
    steps = int(steps)
    reference = bundle.record_map["reference"]
    H = reference.H_spin_free + (reference.H_soc if soc_enabled else 0.0)
    nstate = H.shape[0]
    if checkpoint is None:
        coefficients = _validate_state_v240(initial_coefficients, nstate)
        start_step = 0
        start_time = 0.0
    else:
        if type(checkpoint) is not FrozenSnapshotCheckpointV240:
            raise TypeError("checkpoint must be FrozenSnapshotCheckpointV240.")
        if checkpoint.bundle_fingerprint != bundle.fingerprint:
            raise ValueError("checkpoint and snapshot bundle identities differ.")
        coefficients = np.asarray(checkpoint.coefficients, dtype=complex)
        if coefficients.shape != (nstate,) or not np.all(np.isfinite(coefficients)):
            raise ValueError("checkpoint electronic state is incompatible or non-finite.")
        if abs(float(np.linalg.norm(coefficients)) - 1.0) > 1.0e-10:
            raise ValueError("checkpoint electronic state is not normalized.")
        coefficients = coefficients.copy()
        start_step = int(checkpoint.step)
        start_time = float(checkpoint.time_au)
        if abs(start_time - start_step * dt) > 1.0e-12:
            raise ValueError("checkpoint step/time is inconsistent with the timestep.")
    propagator = expm(-1j * dt * H)
    history = [coefficients.copy()]
    for _ in range(steps):
        coefficients = propagator @ coefficients
        history.append(coefficients.copy())
    history = np.asarray(history, dtype=complex)
    times = start_time + dt * np.arange(steps + 1, dtype=float)
    norms = np.sum(np.abs(history) ** 2, axis=1)
    final_step = start_step + steps
    final_checkpoint = FrozenSnapshotCheckpointV240(
        bundle_fingerprint=bundle.fingerprint,
        step=final_step,
        time_au=float(final_step * dt),
        coefficients=history[-1].copy(),
    )
    return FrozenSnapshotDynamicsV240(
        bundle_fingerprint=bundle.fingerprint,
        evidence_class=bundle.source_kind,
        soc_enabled=soc_enabled,
        times_au=times,
        coefficients=history,
        norms=norms,
        checkpoint=final_checkpoint,
    )


def run_admitted_external_soc_dynamics_v240(
    bundle,
    admission_audit,
    initial_coefficients,
    **kwargs,
):
    if type(admission_audit) is not ExternalSOCAdmissionAuditV240:
        raise TypeError("external dynamics requires the typed v0.24.0 admission audit.")
    if not admission_audit.external_snapshot_admitted:
        raise ValueError("external SOC dynamics requires an admitted external snapshot.")
    if admission_audit.bundle_fingerprint != bundle.fingerprint:
        raise ValueError("admission proof and dynamics bundle fingerprints differ.")
    return preview_frozen_snapshot_dynamics_v240(
        bundle, initial_coefficients, **kwargs
    )
