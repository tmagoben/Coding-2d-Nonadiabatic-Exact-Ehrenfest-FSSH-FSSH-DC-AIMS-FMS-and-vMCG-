from dataclasses import dataclass, field
import numpy as np


@dataclass
class LocalDiabaticTBF:
    """Gaussian TBF carrying an explicit electronic spinor.

    `state` is the adiabatic surface used to guide the classical nuclear center.
    `spinor` is the electronic vector used in the Gaussian basis function and is
    represented in the analytic model's global diabatic basis.

    Keeping these concepts separate is essential in v0.12:
      - guidance state controls qdot/pdot;
      - electronic spinor controls quantum S/H/T matrix elements.
    """
    uid: int
    state: int
    q: np.ndarray
    p: np.ndarray
    A: np.ndarray
    spinor: np.ndarray
    node: object = None
    spawned_targets: set = field(default_factory=set)

    def __post_init__(self):
        self.uid = int(self.uid)
        self.state = int(self.state)
        self.q = np.asarray(self.q, dtype=float)
        self.p = np.asarray(self.p, dtype=float)
        self.A = np.asarray(self.A, dtype=float)
        self.spinor = np.asarray(self.spinor, dtype=complex)

        if self.q.ndim != 1 or self.p.shape != self.q.shape:
            raise ValueError("q and p must be equal-length vectors.")
        if self.A.shape != (len(self.q), len(self.q)):
            raise ValueError("A has incompatible shape.")
        if self.spinor.ndim != 1:
            raise ValueError("spinor must be one-dimensional.")

        n = np.linalg.norm(self.spinor)
        if n <= 0.0:
            raise ValueError("spinor cannot be zero.")
        self.spinor = self.spinor/n

    def copy(self):
        return LocalDiabaticTBF(
            uid=self.uid,
            state=self.state,
            q=self.q.copy(),
            p=self.p.copy(),
            A=self.A.copy(),
            spinor=self.spinor.copy(),
            node=self.node,
            spawned_targets=set(self.spawned_targets),
        )


def from_adiabatic_guided_tbf(tbf, provider):
    point = provider.evaluate(np.asarray(tbf.q, dtype=float))
    spinor = np.asarray(point.frame[:, int(tbf.state)], dtype=complex)
    return LocalDiabaticTBF(
        uid=tbf.uid,
        state=tbf.state,
        q=tbf.q,
        p=tbf.p,
        A=tbf.A,
        spinor=spinor,
        node=getattr(tbf, "node", None),
        spawned_targets=set(getattr(tbf, "spawned_targets", set())),
    )


def reset_to_instantaneous_adiabatic_spinor(tbf, provider):
    point = provider.evaluate(np.asarray(tbf.q, dtype=float))
    tbf.spinor = np.asarray(
        point.frame[:, int(tbf.state)],
        dtype=complex,
    )
    return tbf


def parallel_transport_spinor_full_space(tbf, old_q, new_q, provider):
    r"""Parallel transport through complete old/new adiabatic frames.

    Let
        O = Phi_old^dag Phi_new.

    If c_old are coordinates of the physical spinor in Phi_old, then its coordinates
    in Phi_new are
        c_new = O^dag c_old.

    For a complete electronic frame this leaves the physical vector in the global
    diabatic basis unchanged exactly.  The function is written explicitly to expose
    the overlap-transport algebra used later for ab initio/local-diabatic extensions.
    """
    old_point = provider.evaluate(np.asarray(old_q, dtype=float))
    new_point = provider.evaluate(np.asarray(new_q, dtype=float))

    Phi_old = np.asarray(old_point.frame, dtype=complex)
    Phi_new = np.asarray(new_point.frame, dtype=complex)

    c_old = Phi_old.conj().T @ np.asarray(tbf.spinor, dtype=complex)
    O = Phi_old.conj().T @ Phi_new
    c_new = O.conj().T @ c_old
    transported = Phi_new @ c_new

    transported /= np.linalg.norm(transported)
    tbf.spinor = transported
    return tbf
