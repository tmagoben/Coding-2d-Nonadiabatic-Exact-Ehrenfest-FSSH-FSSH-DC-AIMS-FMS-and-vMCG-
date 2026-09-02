from dataclasses import dataclass
import numpy as np

from .gaussian_nd import analytic_overlap_equal_width


def _validate_spd(M, name):
    M = np.asarray(M, dtype=float)
    if M.ndim != 2 or M.shape[0] != M.shape[1]:
        raise ValueError(f"{name} must be square.")
    if not np.allclose(M, M.T, atol=1e-12):
        raise ValueError(f"{name} must be symmetric.")
    if np.min(np.linalg.eigvalsh(M)) <= 0.0:
        raise ValueError(f"{name} must be positive definite.")
    return M


def overlap_centroid_equal_width(qi, pi, qj, pj, A):
    A = _validate_spd(A, "A")
    qi = np.asarray(qi, float)
    qj = np.asarray(qj, float)
    pi = np.asarray(pi, float)
    pj = np.asarray(pj, float)
    return 0.5*(qi+qj) + 0.5j*np.linalg.solve(A, pj-pi)


def gradient_matrix_element_equal_width(qi, pi, qj, pj, A):
    """Return <g_i|grad g_j> as a D-vector."""
    A = _validate_spd(A, "A")
    S = analytic_overlap_equal_width(qi, pi, qj, pj, A)
    mu = overlap_centroid_equal_width(qi, pi, qj, pj, A)
    factor = -A @ (mu - np.asarray(qj, float)) + 1j*np.asarray(pj, float)
    return factor * S


def kinetic_matrix_element_equal_width(qi, pi, qj, pj, A, mass_matrix):
    """Return <g_i| -1/2 grad^T M^-1 grad |g_j>."""
    A = _validate_spd(A, "A")
    M = _validate_spd(mass_matrix, "mass_matrix")
    B = np.linalg.inv(M)

    S = analytic_overlap_equal_width(qi, pi, qj, pj, A)
    mu = overlap_centroid_equal_width(qi, pi, qj, pj, A)

    ui = -A @ (mu - np.asarray(qi, float)) - 1j*np.asarray(pi, float)
    uj = -A @ (mu - np.asarray(qj, float)) + 1j*np.asarray(pj, float)

    return 0.5 * S * (
        ui @ B @ uj
        + 0.5*np.trace(B @ A)
    )


def basis_time_matrix_element_equal_width(
    qi, pi, qj, pj, A, qdot_j, pdot_j
):
    """Return <g_i|dot g_j> for fixed-width zero-chirp Gaussians."""
    A = _validate_spd(A, "A")
    S = analytic_overlap_equal_width(qi, pi, qj, pj, A)
    mu = overlap_centroid_equal_width(qi, pi, qj, pj, A)

    y = mu - np.asarray(qj, float)
    factor = (
        (A @ y - 1j*np.asarray(pj, float)) @ np.asarray(qdot_j, float)
        + 1j*y @ np.asarray(pdot_j, float)
    )
    return S * factor


def local_d2_matrix(nac_q, mass_matrix):
    """Electronic matrix D2_ab = sum_c,alpha,beta d_ac,a M^-1_ab d_cb,b."""
    d = np.asarray(nac_q, dtype=float)
    M = _validate_spd(mass_matrix, "mass_matrix")
    B = np.linalg.inv(M)

    ns, ns2, ndim = d.shape
    if ns != ns2 or M.shape != (ndim, ndim):
        raise ValueError("NAC tensor and mass matrix have incompatible shapes.")

    return np.einsum("ack,kl,cbl->ab", d, B, d)


@dataclass
class LocalAdiabaticTBF:
    state: int
    q: np.ndarray
    p: np.ndarray
    A: np.ndarray

    def __post_init__(self):
        self.q = np.asarray(self.q, dtype=float)
        self.p = np.asarray(self.p, dtype=float)
        self.A = _validate_spd(self.A, "A")
        if self.q.shape != self.p.shape or self.A.shape != (len(self.q), len(self.q)):
            raise ValueError("TBF q, p, and A shapes are inconsistent.")

    def copy(self):
        return LocalAdiabaticTBF(
            int(self.state),
            self.q.copy(),
            self.p.copy(),
            self.A.copy(),
        )


def local_pair_hamiltonian_element(tbf_i, tbf_j, provider):
    """Local constant-electronic-quantity matrix element at pair midpoint.

    Includes:
      - exact equal-width Gaussian kinetic integral,
      - centroid adiabatic energy,
      - first-order NAC kinetic coupling,
      - locally constant d^2 term.

    Neglects:
      - derivative/divergence of the NAC field over the overlap region.
    """
    if not np.allclose(tbf_i.A, tbf_j.A, atol=1e-12):
        raise ValueError("v0.5 local pair matrix elements require equal width matrices.")

    qbar = 0.5*(tbf_i.q + tbf_j.q)
    point = provider.evaluate(qbar)
    M = point.mass_matrix_q_au
    B = np.linalg.inv(M)

    S_nuclear = analytic_overlap_equal_width(
        tbf_i.q, tbf_i.p, tbf_j.q, tbf_j.p, tbf_i.A
    )
    G = gradient_matrix_element_equal_width(
        tbf_i.q, tbf_i.p, tbf_j.q, tbf_j.p, tbf_i.A
    )
    Tij = kinetic_matrix_element_equal_width(
        tbf_i.q, tbf_i.p, tbf_j.q, tbf_j.p, tbf_i.A, M
    )

    a = tbf_i.state
    b = tbf_j.state

    value = 0.0 + 0.0j

    if a == b:
        value += Tij + point.energies[a] * S_nuclear

    d_ab = point.nac_q[a, b]
    value += -d_ab @ B @ G

    D2 = local_d2_matrix(point.nac_q, M)
    value += -0.5 * D2[a, b] * S_nuclear

    return value


def local_overlap_element(tbf_i, tbf_j):
    if tbf_i.state != tbf_j.state:
        return 0.0 + 0.0j
    if not np.allclose(tbf_i.A, tbf_j.A, atol=1e-12):
        raise ValueError("v0.5 local overlap requires equal widths.")
    return analytic_overlap_equal_width(
        tbf_i.q, tbf_i.p, tbf_j.q, tbf_j.p, tbf_i.A
    )


def tbf_guidance(tbf, provider):
    point = provider.evaluate(tbf.q)
    M = point.mass_matrix_q_au
    qdot = np.linalg.solve(M, tbf.p)
    pdot = -point.gradients_q[tbf.state]
    return qdot, pdot


def local_basis_time_element(tbf_i, tbf_j, provider):
    if tbf_i.state != tbf_j.state:
        return 0.0 + 0.0j
    if not np.allclose(tbf_i.A, tbf_j.A, atol=1e-12):
        raise ValueError("v0.5 local basis-time matrix requires equal widths.")

    qdot_j, pdot_j = tbf_guidance(tbf_j, provider)

    return basis_time_matrix_element_equal_width(
        tbf_i.q,
        tbf_i.p,
        tbf_j.q,
        tbf_j.p,
        tbf_j.A,
        qdot_j,
        pdot_j,
    )


def local_matrices(basis, provider):
    n = len(basis)
    S = np.zeros((n, n), dtype=complex)
    H = np.zeros((n, n), dtype=complex)
    T = np.zeros((n, n), dtype=complex)

    for i in range(n):
        for j in range(n):
            S[i, j] = local_overlap_element(basis[i], basis[j])
            H[i, j] = local_pair_hamiltonian_element(basis[i], basis[j], provider)
            T[i, j] = local_basis_time_element(basis[i], basis[j], provider)

    return S, H, T
