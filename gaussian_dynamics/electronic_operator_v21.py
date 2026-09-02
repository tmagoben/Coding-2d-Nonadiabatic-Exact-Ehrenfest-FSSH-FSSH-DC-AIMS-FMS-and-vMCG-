from dataclasses import dataclass, field
import numpy as np

from .matrix_invariants_v213 import (
    antihermiticity_residual_v213,
    hermiticity_residual_v213,
    isometry_residual_v213,
    require_residual_v213,
    symmetry_residual_v213,
)


def _strict_real_array_v213(value, name):
    raw = np.asarray(value)
    if np.iscomplexobj(raw) and np.any(np.imag(raw) != 0.0):
        raise ValueError(f"{name} must be real; complex data would be discarded.")
    return np.asarray(np.real(raw), dtype=float)


@dataclass
class ElectronicOperatorPointV21:
    q: np.ndarray
    H: np.ndarray
    dH_dq: np.ndarray
    connection_q: np.ndarray
    mass_matrix_q_au: np.ndarray
    metadata: dict = field(default_factory=dict)

    def structural_residuals_v213(self):
        H = np.asarray(self.H, dtype=complex)
        dH = np.asarray(self.dH_dq, dtype=complex)
        D = np.asarray(self.connection_q, dtype=complex)
        M = _strict_real_array_v213(self.mass_matrix_q_au, "mass_matrix_q_au")
        return {
            "H_hermiticity": hermiticity_residual_v213(H),
            "dH_hermiticity": tuple(
                hermiticity_residual_v213(dH[a]) for a in range(len(dH))
            ),
            "connection_antihermiticity": tuple(
                antihermiticity_residual_v213(D[a]) for a in range(len(D))
            ),
            "mass_symmetry": symmetry_residual_v213(M),
        }

    def validate(self, atol=1e-12):
        atol = float(atol)
        if not np.isfinite(atol) or atol < 0.0:
            raise ValueError("operator structural tolerance must be finite and nonnegative.")
        q = _strict_real_array_v213(self.q, "q")
        H = np.asarray(self.H, dtype=complex)
        dH = np.asarray(self.dH_dq, dtype=complex)
        D = np.asarray(self.connection_q, dtype=complex)
        M = _strict_real_array_v213(self.mass_matrix_q_au, "mass_matrix_q_au")
        if q.ndim != 1:
            raise ValueError("q must be a one-dimensional coordinate vector.")
        if H.ndim != 2 or H.shape[0] != H.shape[1]:
            raise ValueError("H must be a square matrix.")
        if H.shape[0] == 0:
            raise ValueError("H must contain at least one electronic state.")
        ns = H.shape[0]
        nq = len(q)
        if dH.shape != (nq, ns, ns):
            raise ValueError("dH_dq must have shape (nq,nstate,nstate).")
        if D.shape != (nq, ns, ns):
            raise ValueError("connection_q must have shape (nq,nstate,nstate).")
        if M.shape != (nq, nq):
            raise ValueError("mass_matrix_q_au must have shape (nq,nq).")
        for arr, name in ((q, "q"), (H, "H"), (dH, "dH_dq"), (D, "connection_q"), (M, "mass_matrix_q_au")):
            if not np.all(np.isfinite(arr)):
                raise ValueError(f"{name} contains non-finite data.")
        require_residual_v213(
            "Electronic Hamiltonian Hermiticity",
            hermiticity_residual_v213(H),
            atol,
        )
        for a in range(nq):
            require_residual_v213(
                f"Hamiltonian-derivative Hermiticity at coordinate {a}",
                hermiticity_residual_v213(dH[a]),
                atol,
            )
            require_residual_v213(
                f"derivative-connection anti-Hermiticity at coordinate {a}",
                antihermiticity_residual_v213(D[a]),
                atol,
            )
        require_residual_v213(
            "generalized-mass symmetry",
            symmetry_residual_v213(M),
            atol,
        )
        if np.min(np.linalg.eigvalsh(M)) <= 0.0:
            raise ValueError("mass_matrix_q_au must be positive definite.")
        self.q, self.H, self.dH_dq, self.connection_q, self.mass_matrix_q_au = q, H, dH, D, M
        return self

    @property
    def nstate(self):
        return int(self.H.shape[0])

    @property
    def nq(self):
        return int(len(self.q))

    @property
    def hamiltonian_derivative_operator_q(self):
        """Unambiguous alias for the physical operator matrices ``dH_dq``.

        These matrices are ``<Phi|partial_a H_operator|Phi>``.  They transform
        covariantly and are not the naive coordinate derivative of a matrix written in
        an arbitrarily moving electronic frame.
        """
        return self.dH_dq

    def force_expectation(self, electronic_vector):
        c = np.asarray(electronic_vector, dtype=complex)
        if c.shape != (self.nstate,):
            raise ValueError("electronic_vector has incompatible shape.")
        if not np.all(np.isfinite(c)):
            raise ValueError("electronic_vector contains non-finite data.")
        norm = float(np.real(np.vdot(c, c)))
        if norm <= 0.0:
            raise ValueError("electronic_vector must have nonzero norm.")
        c = c / np.sqrt(norm)
        return np.asarray([-np.real(np.vdot(c, self.dH_dq[a] @ c)) for a in range(self.nq)], dtype=float)


@dataclass
class ElectronicOperatorSnapshotV21:
    point: ElectronicOperatorPointV21
    state_vectors: np.ndarray | None = None
    wavefunction_snapshot: object | None = None
    parent_snapshot: object | None = None
    frame_from_parent: np.ndarray | None = None
    metadata: dict = field(default_factory=dict)

    def validate(self, atol=1e-12, isometry_atol=1e-10):
        self.point = self.point.validate(atol=atol)
        ns = self.point.nstate
        if self.state_vectors is not None:
            V = np.asarray(self.state_vectors, dtype=complex)
            if V.ndim != 2 or V.shape[1] != ns:
                raise ValueError("state_vectors must have shape (representation_dimension,nstate).")
            require_residual_v213(
                "state-vector isometry",
                isometry_residual_v213(V),
                isometry_atol,
            )
            self.state_vectors = V
        if self.frame_from_parent is not None:
            G = np.asarray(self.frame_from_parent, dtype=complex)
            if G.shape != (ns, ns):
                raise ValueError("frame_from_parent has incompatible shape.")
            require_residual_v213(
                "frame_from_parent unitarity",
                isometry_residual_v213(G),
                isometry_atol,
            )
            self.frame_from_parent = G
        return self


def adiabatic_point_to_operator_v21(point):
    point = point.validate()
    E = np.asarray(point.energies, dtype=float)
    grad = np.asarray(point.gradients_q, dtype=float)
    nac = np.asarray(point.nac_q, dtype=float)
    ns = len(E)
    nq = len(point.q)
    H = np.diag(E).astype(complex)
    dH = np.zeros((nq, ns, ns), dtype=complex)
    D = np.zeros((nq, ns, ns), dtype=complex)
    for a in range(nq):
        D[a] = nac[:, :, a]
        dH[a][np.diag_indices(ns)] = grad[:, a]
        for i in range(ns):
            for j in range(ns):
                if i != j:
                    dH[a, i, j] = (E[j] - E[i]) * nac[i, j, a]
    return ElectronicOperatorPointV21(
        q=np.asarray(point.q, float).copy(), H=H, dH_dq=dH, connection_q=D,
        mass_matrix_q_au=np.asarray(point.mass_matrix_q_au, float).copy(),
        metadata={**dict(point.metadata), "v21_operator_conversion": "adiabatic_real_to_full_complex"},
    ).validate()


class ElectronicOperatorProviderAdapterV21:
    def __init__(self, base_provider):
        self.base_provider = base_provider

    def evaluate_snapshot(self, q):
        base = self.base_provider.evaluate_snapshot(q)
        point = adiabatic_point_to_operator_v21(base.point)
        ns = point.nstate
        return ElectronicOperatorSnapshotV21(
            point=point,
            state_vectors=None if base.state_vectors is None else np.asarray(base.state_vectors, complex).copy(),
            wavefunction_snapshot=base.wavefunction_snapshot,
            parent_snapshot=base,
            frame_from_parent=np.eye(ns, dtype=complex),
            metadata={"adapter": "ElectronicOperatorProviderAdapterV21", "source_node_id": base.node_id},
        ).validate()

    def evaluate(self, q):
        return self.evaluate_snapshot(q).point

    def snapshot_overlap(self, left, right):
        if left.state_vectors is not None and right.state_vectors is not None:
            return left.state_vectors.conj().T @ right.state_vectors
        if left.parent_snapshot is not None and right.parent_snapshot is not None and hasattr(self.base_provider, "snapshot_overlap"):
            O = self.base_provider.snapshot_overlap(left.parent_snapshot, right.parent_snapshot)
            Gl = np.eye(left.point.nstate, dtype=complex) if left.frame_from_parent is None else left.frame_from_parent
            Gr = np.eye(right.point.nstate, dtype=complex) if right.frame_from_parent is None else right.frame_from_parent
            return Gl.conj().T @ O @ Gr
        raise ValueError("No cross-geometry overlap path is available.")

    def diagnostics_dict(self):
        return self.base_provider.diagnostics_dict() if hasattr(self.base_provider, "diagnostics_dict") else {}
