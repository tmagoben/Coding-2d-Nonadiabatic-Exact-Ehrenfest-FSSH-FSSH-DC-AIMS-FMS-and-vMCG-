from dataclasses import dataclass
import numpy as np
from .electronic_operator_v21 import ElectronicOperatorPointV21, ElectronicOperatorSnapshotV21
from .matrix_invariants_v213 import isometry_residual_v213, require_residual_v213


def random_unitary_v21(n, seed=210021):
    rng = np.random.default_rng(int(seed))
    X = rng.normal(size=(int(n), int(n))) + 1j * rng.normal(size=(int(n), int(n)))
    Q, R = np.linalg.qr(X)
    d = np.diag(R)
    phase = np.where(np.abs(d) > 1e-15, d / np.abs(d), 1.0)
    return Q * phase.conj()[None, :]


@dataclass(frozen=True)
class PhaseMixingGaugeV21:
    U0: np.ndarray
    phase_gradient: np.ndarray
    phase_offset: np.ndarray

    def __post_init__(self):
        U = np.asarray(self.U0, dtype=complex)
        B = np.asarray(self.phase_gradient, dtype=float)
        t0 = np.asarray(self.phase_offset, dtype=float)
        if U.ndim != 2 or U.shape[0] != U.shape[1]:
            raise ValueError("U0 must be square.")
        ns = U.shape[0]
        require_residual_v213(
            "U0 unitarity", isometry_residual_v213(U), 1.0e-10
        )
        if B.ndim != 2 or B.shape[0] != ns:
            raise ValueError("phase_gradient must have shape (nstate,nq).")
        if t0.shape != (ns,):
            raise ValueError("phase_offset has incompatible shape.")
        object.__setattr__(self, "U0", U)
        object.__setattr__(self, "phase_gradient", B)
        object.__setattr__(self, "phase_offset", t0)

    @property
    def nstate(self): return self.U0.shape[0]
    @property
    def nq(self): return self.phase_gradient.shape[1]

    def matrix(self, q):
        q = np.asarray(q, float)
        theta = self.phase_offset + self.phase_gradient @ q
        return np.diag(np.exp(1j * theta)) @ self.U0

    def derivatives(self, q):
        q = np.asarray(q, float)
        theta = self.phase_offset + self.phase_gradient @ q
        phase = np.exp(1j * theta)
        out = np.empty((self.nq, self.nstate, self.nstate), dtype=complex)
        for a in range(self.nq):
            out[a] = np.diag(1j * self.phase_gradient[:, a] * phase) @ self.U0
        return out

    def velocity_derivative(self, q, qdot):
        return np.tensordot(np.asarray(qdot, float), self.derivatives(q), axes=(0, 0))


def transform_operator_point_v21(point, G, dG_dq):
    point = point.validate()
    G = np.asarray(G, complex)
    dG = np.asarray(dG_dq, complex)
    H = G.conj().T @ point.H @ G
    dH = np.asarray([G.conj().T @ point.dH_dq[a] @ G for a in range(point.nq)])
    D = np.asarray([G.conj().T @ point.connection_q[a] @ G + G.conj().T @ dG[a] for a in range(point.nq)])
    return ElectronicOperatorPointV21(point.q.copy(), H, dH, D, point.mass_matrix_q_au.copy(), {**dict(point.metadata), "v21_complex_gauge": True}).validate()


class GaugeTransformedOperatorProviderV21:
    def __init__(self, base_provider, gauge):
        self.base_provider = base_provider
        self.gauge = gauge

    def evaluate_snapshot(self, q):
        base = self.base_provider.evaluate_snapshot(q)
        G = self.gauge.matrix(q)
        point = transform_operator_point_v21(base.point, G, self.gauge.derivatives(q))
        V = None if base.state_vectors is None else base.state_vectors @ G
        return ElectronicOperatorSnapshotV21(point, V, base.wavefunction_snapshot, base, G, {**dict(base.metadata), "gauge_provider": "GaugeTransformedOperatorProviderV21"}).validate()

    def evaluate(self, q): return self.evaluate_snapshot(q).point

    def snapshot_overlap(self, left, right):
        if left.state_vectors is not None and right.state_vectors is not None:
            return left.state_vectors.conj().T @ right.state_vectors
        O = self.base_provider.snapshot_overlap(left.parent_snapshot, right.parent_snapshot)
        return left.frame_from_parent.conj().T @ O @ right.frame_from_parent

    def diagnostics_dict(self):
        return self.base_provider.diagnostics_dict() if hasattr(self.base_provider, "diagnostics_dict") else {}
