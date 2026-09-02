import numpy as np

from .gaussian import frozen_gaussian, kinetic_on_gaussian
from .grids import inner_product


def basis_functions(x, q, p, alpha):
    """Return matrix G[x_index, gaussian_index]."""
    return np.column_stack([
        frozen_gaussian(x, qj, pj, alpha)
        for qj, pj in zip(q, p)
    ])


def basis_time_derivatives(x, q, p, alpha, qdot, pdot):
    """Time derivative of each fixed-width moving Gaussian."""
    columns = []

    for qj, pj, qdj, pdj in zip(q, p, qdot, pdot):
        g = frozen_gaussian(x, qj, pj, alpha)
        y = x - qj

        dg_dq = (alpha * y - 1j * pj) * g
        dg_dp = 1j * y * g

        columns.append(dg_dq * qdj + dg_dp * pdj)

    return np.column_stack(columns)


def moving_basis_matrices(x, dx, q, p, alpha, mass, potential, gradient):
    """Construct S, H, tau for a classically guided frozen-Gaussian basis."""
    q = np.asarray(q, dtype=float)
    p = np.asarray(p, dtype=float)

    G = basis_functions(x, q, p, alpha)

    qdot = p / mass
    pdot = -np.asarray([gradient(qj) for qj in q], dtype=float)

    Gdot = basis_time_derivatives(x, q, p, alpha, qdot, pdot)

    n = len(q)
    S = np.zeros((n, n), dtype=complex)
    H = np.zeros((n, n), dtype=complex)
    tau = np.zeros((n, n), dtype=complex)

    V = np.asarray(potential(x), dtype=float)

    H_on_basis = []
    for j in range(n):
        Tg = kinetic_on_gaussian(x, q[j], p[j], alpha, mass)
        H_on_basis.append(Tg + V * G[:, j])

    H_on_basis = np.column_stack(H_on_basis)

    for i in range(n):
        for j in range(n):
            S[i, j] = inner_product(G[:, i], G[:, j], dx)
            H[i, j] = inner_product(G[:, i], H_on_basis[:, j], dx)
            tau[i, j] = inner_product(G[:, i], Gdot[:, j], dx)

    return S, H, tau, qdot, pdot


def coefficient_rhs(C, S, H, tau):
    """Solve i S Cdot = (H - i tau) C without explicitly inverting S."""
    rhs = -1j * (H @ C) - tau @ C
    return np.linalg.solve(S, rhs)


def wavefunction_from_basis(x, q, p, alpha, C):
    G = basis_functions(x, q, p, alpha)
    return G @ np.asarray(C, dtype=complex)


def basis_norm(C, S):
    return float(np.real(np.vdot(C, S @ C)))


def _pack(C, q, p):
    return np.concatenate([
        np.asarray(C, dtype=complex),
        np.asarray(q, dtype=complex),
        np.asarray(p, dtype=complex),
    ])


def _unpack(z, n):
    C = z[:n]
    q = z[n:2*n].real
    p = z[2*n:3*n].real
    return C, q, p


def run_moving_gaussian_basis(
    x,
    dx,
    q0,
    p0,
    C0,
    alpha,
    mass,
    potential,
    gradient,
    dt=0.001,
    steps=1000,
    store_every=10,
):
    """Propagate a small moving frozen-Gaussian basis.

    Gaussian centers follow classical equations.
    Coefficients obey the projected TDSE in the nonorthogonal moving basis.
    """
    q0 = np.asarray(q0, dtype=float)
    p0 = np.asarray(p0, dtype=float)
    C0 = np.asarray(C0, dtype=complex)
    n = len(q0)

    if not (len(p0) == len(C0) == n):
        raise ValueError("q0, p0, and C0 must have the same length.")

    S0, _, _, _, _ = moving_basis_matrices(
        x, dx, q0, p0, alpha, mass, potential, gradient
    )
    nrm = np.real(np.vdot(C0, S0 @ C0))
    if nrm <= 0.0:
        raise ValueError("Initial moving-basis wavefunction has zero norm.")
    C0 = C0 / np.sqrt(nrm)

    z = _pack(C0, q0, p0)

    def rhs(zlocal):
        C, q, p = _unpack(zlocal, n)
        S, H, tau, qdot, pdot = moving_basis_matrices(
            x, dx, q, p, alpha, mass, potential, gradient
        )
        Cdot = coefficient_rhs(C, S, H, tau)
        return _pack(Cdot, qdot, pdot)

    times = []
    qs = []
    ps = []
    coeffs = []
    norms = []
    condition_numbers = []
    states = []

    def record(step):
        C, q, p = _unpack(z, n)
        S, _, _, _, _ = moving_basis_matrices(
            x, dx, q, p, alpha, mass, potential, gradient
        )

        times.append(step * dt)
        qs.append(q.copy())
        ps.append(p.copy())
        coeffs.append(C.copy())
        norms.append(basis_norm(C, S))
        condition_numbers.append(np.linalg.cond(S))
        states.append(wavefunction_from_basis(x, q, p, alpha, C))

    record(0)

    for step in range(1, steps + 1):
        k1 = rhs(z)
        k2 = rhs(z + 0.5 * dt * k1)
        k3 = rhs(z + 0.5 * dt * k2)
        k4 = rhs(z + dt * k3)

        z = z + dt * (k1 + 2*k2 + 2*k3 + k4) / 6.0

        if step % store_every == 0:
            record(step)

    return {
        "time": np.asarray(times),
        "q": np.asarray(qs),
        "p": np.asarray(ps),
        "C": np.asarray(coeffs),
        "norm": np.asarray(norms),
        "condition_number": np.asarray(condition_numbers),
        "psi": np.asarray(states),
    }
