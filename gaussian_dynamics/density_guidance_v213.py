"""Degeneracy-safe, representation-covariant nuclear guidance.

No eigenvector is selected when a local electronic coefficient block is empty.  Each
Gaussian instead carries an auxiliary electronic density matrix.  A populated block
refreshes that density from its coefficients; a weak block retains and transports its
last valid density.  A genuinely unseeded zero block receives zero force until a guide
density is inherited or physical amplitude arrives.
"""

from dataclasses import dataclass, asdict
import numpy as np

from .gauge_graph import nearest_unitary
from .matrix_invariants_v213 import (
    hermiticity_residual_v213,
    require_residual_v213,
)


@dataclass(frozen=True)
class DensityMatrixGuidanceSettingsV213:
    minimum_local_amplitude: float = 1.0e-12
    low_amplitude_policy: str = "retained_density"  # retained_density | zero_force
    density_tolerance: float = 1.0e-10

    def validate(self):
        if not np.isfinite(self.minimum_local_amplitude) or self.minimum_local_amplitude < 0.0:
            raise ValueError("minimum_local_amplitude cannot be negative.")
        if self.low_amplitude_policy not in {"retained_density", "zero_force"}:
            raise ValueError("invalid low_amplitude_policy.")
        if not np.isfinite(self.density_tolerance) or self.density_tolerance <= 0.0:
            raise ValueError("density_tolerance must be positive.")
        return self


def normalized_density_from_vector_v213(vector):
    c = np.asarray(vector, dtype=complex)
    if c.ndim != 1:
        raise ValueError("electronic vector must be one-dimensional.")
    if not np.all(np.isfinite(c)):
        raise ValueError("electronic vector contains non-finite data.")
    norm = float(np.real(np.vdot(c, c)))
    if norm <= 0.0:
        raise ValueError("electronic vector must have nonzero norm.")
    return np.outer(c, c.conj()) / norm


def validate_guide_density_v213(density, nstate, tolerance=1.0e-10, *, allow_zero=True):
    rho = np.asarray(density, dtype=complex)
    if rho.shape != (int(nstate), int(nstate)):
        raise ValueError("guide density has incompatible shape.")
    if not np.all(np.isfinite(rho)):
        raise ValueError("guide density contains non-finite data.")
    require_residual_v213(
        "guide-density Hermiticity",
        hermiticity_residual_v213(rho),
        tolerance,
    )
    rho = 0.5 * (rho + rho.conj().T)
    trace = float(np.real(np.trace(rho)))
    if allow_zero and abs(trace) <= tolerance:
        if np.linalg.norm(rho, ord="fro") > tolerance:
            raise ValueError("a zero-trace guide density must be the zero matrix.")
        return np.zeros_like(rho)
    if abs(trace - 1.0) > tolerance:
        raise ValueError("a nonzero guide density must have unit trace.")
    minimum = float(np.min(np.linalg.eigvalsh(rho)))
    if minimum < -tolerance:
        raise ValueError("guide density must be positive semidefinite.")
    return rho / trace


def density_force_v213(density, derivative_operators):
    rho = np.asarray(density, dtype=complex)
    K = np.asarray(derivative_operators, dtype=complex)
    if K.ndim != 3 or K.shape[1:] != rho.shape:
        raise ValueError("Hamiltonian-derivative operators are incompatible with density.")
    return np.asarray(
        [-float(np.real(np.trace(rho @ K[a]))) for a in range(K.shape[0])],
        dtype=float,
    )


class BlockDensityMatrixGuidanceV213:
    """Stateful guide-density manager keyed by stable Gaussian uid."""

    def __init__(self, settings=DensityMatrixGuidanceSettingsV213()):
        self.settings = settings.validate()
        self._densities = {}
        self._snapshots = {}
        self.coefficient_refreshes = 0
        self.density_transports = 0
        self.retained_density_uses = 0
        self.zero_force_uses = 0
        self.parent_inheritances = 0
        self.explicit_seeds = 0

    def _transport_to_snapshot(self, uid, snapshot, provider):
        uid = int(uid)
        if uid not in self._densities or uid not in self._snapshots:
            return
        previous = self._snapshots[uid]
        overlap = np.asarray(provider.snapshot_overlap(previous, snapshot), dtype=complex)
        if overlap.shape != (snapshot.point.nstate, snapshot.point.nstate):
            raise ValueError("snapshot overlap has incompatible electronic dimension.")
        U = nearest_unitary(overlap)
        rho = U.conj().T @ self._densities[uid] @ U
        self._densities[uid] = validate_guide_density_v213(
            rho,
            snapshot.point.nstate,
            self.settings.density_tolerance,
        )
        self.density_transports += 1

    def forces_and_masses(self, basis, coefficients, provider, nstate):
        basis = list(basis)
        C = np.asarray(coefficients, dtype=complex)
        s = int(nstate)
        if s < 1 or s != nstate:
            raise ValueError("nstate must be a positive integer.")
        if C.shape != (len(basis) * s,):
            raise ValueError("coefficient dimension is inconsistent with basis*nstate.")
        if not np.all(np.isfinite(C)):
            raise ValueError("electronic coefficients contain non-finite data.")
        uids = [int(item.uid) for item in basis]
        if len(set(uids)) != len(uids):
            raise ValueError("Gaussian uids must be unique for density guidance.")

        forces = []
        masses = []
        local_norms = []
        for i, b in enumerate(basis):
            uid = int(b.uid)
            snapshot = provider.evaluate_snapshot(b.q)
            if snapshot.point.nstate != s:
                raise ValueError("electronic model-space dimension changed during guidance.")
            self._transport_to_snapshot(uid, snapshot, provider)

            c = C[s * i : s * (i + 1)]
            amplitude = float(np.real(np.vdot(c, c)))
            local_norms.append(amplitude)
            if amplitude > self.settings.minimum_local_amplitude:
                rho = normalized_density_from_vector_v213(c)
                self._densities[uid] = rho
                self.coefficient_refreshes += 1
                force_density = rho
            elif (
                self.settings.low_amplitude_policy == "retained_density"
                and uid in self._densities
            ):
                force_density = self._densities[uid]
                self.retained_density_uses += 1
            else:
                force_density = np.zeros((s, s), dtype=complex)
                self.zero_force_uses += 1

            self._snapshots[uid] = snapshot
            forces.append(
                density_force_v213(
                    force_density,
                    snapshot.point.hamiltonian_derivative_operator_q,
                )
            )
            masses.append(np.asarray(snapshot.point.mass_matrix_q_au, dtype=float))

        return np.asarray(forces, dtype=float), tuple(masses), np.asarray(local_norms)

    def on_insert(
        self,
        tbf,
        provider,
        *,
        parent_uid=None,
        guide_density=None,
    ):
        uid = int(tbf.uid)
        if uid in self._densities or uid in self._snapshots:
            raise ValueError(f"guidance state already exists for Gaussian uid {uid}.")
        if parent_uid is not None and guide_density is not None:
            raise ValueError("choose parent_uid or guide_density, not both.")
        snapshot = provider.evaluate_snapshot(tbf.q)
        s = snapshot.point.nstate
        if guide_density is not None:
            self._densities[uid] = validate_guide_density_v213(
                guide_density, s, self.settings.density_tolerance
            )
            self.explicit_seeds += 1
        elif parent_uid is not None:
            parent_uid = int(parent_uid)
            if parent_uid not in self._densities or parent_uid not in self._snapshots:
                raise ValueError("parent Gaussian has no guide density to inherit.")
            overlap = np.asarray(
                provider.snapshot_overlap(self._snapshots[parent_uid], snapshot),
                dtype=complex,
            )
            if overlap.shape != (s, s):
                raise ValueError("parent-child overlap has incompatible dimension.")
            U = nearest_unitary(overlap)
            inherited = U.conj().T @ self._densities[parent_uid] @ U
            self._densities[uid] = validate_guide_density_v213(
                inherited, s, self.settings.density_tolerance
            )
            self.parent_inheritances += 1
        self._snapshots[uid] = snapshot

    def on_prune(self, uid):
        uid = int(uid)
        self._densities.pop(uid, None)
        self._snapshots.pop(uid, None)

    def density(self, uid):
        uid = int(uid)
        if uid not in self._densities:
            return None
        return self._densities[uid].copy()

    def diagnostics_dict(self):
        return {
            "settings": asdict(self.settings),
            "tracked_densities": int(len(self._densities)),
            "tracked_snapshots": int(len(self._snapshots)),
            "coefficient_refreshes": int(self.coefficient_refreshes),
            "density_transports": int(self.density_transports),
            "retained_density_uses": int(self.retained_density_uses),
            "zero_force_uses": int(self.zero_force_uses),
            "parent_inheritances": int(self.parent_inheritances),
            "explicit_seeds": int(self.explicit_seeds),
        }

    def checkpoint_state(self):
        """Return an opaque copy used to isolate predictor/corrector trial state."""
        return {
            "densities": {uid: rho.copy() for uid, rho in self._densities.items()},
            "snapshots": dict(self._snapshots),
            "counters": (
                self.coefficient_refreshes,
                self.density_transports,
                self.retained_density_uses,
                self.zero_force_uses,
                self.parent_inheritances,
                self.explicit_seeds,
            ),
        }

    def restore_state(self, checkpoint):
        if not isinstance(checkpoint, dict) or set(checkpoint) != {
            "densities", "snapshots", "counters"
        }:
            raise ValueError("invalid density-guidance checkpoint.")
        self._densities = {
            int(uid): np.asarray(rho, dtype=complex).copy()
            for uid, rho in checkpoint["densities"].items()
        }
        self._snapshots = dict(checkpoint["snapshots"])
        (
            self.coefficient_refreshes,
            self.density_transports,
            self.retained_density_uses,
            self.zero_force_uses,
            self.parent_inheritances,
            self.explicit_seeds,
        ) = checkpoint["counters"]
