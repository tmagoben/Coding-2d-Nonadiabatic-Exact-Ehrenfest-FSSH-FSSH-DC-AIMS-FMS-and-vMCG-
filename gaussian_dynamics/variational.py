import numpy as np

from .gaussian import frozen_gaussian
from .grids import inner_product


PARAMETERS_PER_GAUSSIAN = 6


def pack_parameters(coefficients, q, p, alpha, chirp=None):
    """Pack complex coefficients and Gaussian parameters into a real vector.

    Per Gaussian:
        Re(C), Im(C), q, p, log(alpha), chirp
    """
    coefficients = np.asarray(coefficients, dtype=complex)
    q = np.asarray(q, dtype=float)
    p = np.asarray(p, dtype=float)
    alpha = np.asarray(alpha, dtype=float)

    if chirp is None:
        chirp = np.zeros_like(alpha)
    chirp = np.asarray(chirp, dtype=float)

    n = len(coefficients)
    if not (len(q) == len(p) == len(alpha) == len(chirp) == n):
        raise ValueError("All Gaussian parameter arrays must have the same length.")
    if np.any(alpha <= 0.0):
        raise ValueError("All alpha values must be positive.")

    theta = np.zeros(PARAMETERS_PER_GAUSSIAN * n, dtype=float)

    for j in range(n):
        base = PARAMETERS_PER_GAUSSIAN * j
        theta[base + 0] = coefficients[j].real
        theta[base + 1] = coefficients[j].imag
        theta[base + 2] = q[j]
        theta[base + 3] = p[j]
        theta[base + 4] = np.log(alpha[j])
        theta[base + 5] = chirp[j]

    return theta


def unpack_parameters(theta):
    theta = np.asarray(theta, dtype=float)

    if len(theta) % PARAMETERS_PER_GAUSSIAN != 0:
        raise ValueError("theta has an invalid length.")

    n = len(theta) // PARAMETERS_PER_GAUSSIAN

    C = np.zeros(n, dtype=complex)
    q = np.zeros(n)
    p = np.zeros(n)
    alpha = np.zeros(n)
    chirp = np.zeros(n)

    for j in range(n):
        base = PARAMETERS_PER_GAUSSIAN * j
        C[j] = theta[base + 0] + 1j * theta[base + 1]
        q[j] = theta[base + 2]
        p[j] = theta[base + 3]
        alpha[j] = np.exp(theta[base + 4])
        chirp[j] = theta[base + 5]

    return C, q, p, alpha, chirp


def variational_wavefunction(x, theta):
    C, q, p, alpha, chirp = unpack_parameters(theta)

    psi = np.zeros(len(x), dtype=complex)
    for Cj, qj, pj, aj, kj in zip(C, q, p, alpha, chirp):
        psi += Cj * frozen_gaussian(x, qj, pj, aj, chirp=kj)

    return psi


def hamiltonian_action(psi, x, dx, mass, potential):
    """Apply T+V to an arbitrary grid wavefunction using FFT differentiation."""
    k = 2.0 * np.pi * np.fft.fftfreq(len(x), d=dx)
    second = np.fft.ifft(-(k**2) * np.fft.fft(psi))
    Tpsi = -second / (2.0 * mass)
    return Tpsi + np.asarray(potential(x), dtype=float) * psi


def tangent_vectors(x, theta, relative_step=1.0e-6):
    """Numerical central-difference tangent vectors dPsi/dtheta_mu."""
    theta = np.asarray(theta, dtype=float)
    D = np.zeros((len(x), len(theta)), dtype=complex)

    for mu in range(len(theta)):
        h = relative_step * max(1.0, abs(theta[mu]))

        tp = theta.copy()
        tm = theta.copy()
        tp[mu] += h
        tm[mu] -= h

        D[:, mu] = (
            variational_wavefunction(x, tp)
            - variational_wavefunction(x, tm)
        ) / (2.0 * h)

    return D


def tdvp_velocity(
    x,
    dx,
    theta,
    mass,
    potential,
    relative_step=1.0e-6,
    rcond=1.0e-10,
):
    """McLachlan TDVP velocity for real variational parameters.

    G_{mu,nu} = Re <D_mu|D_nu>
    b_mu       = Im <D_mu|H Psi>
    G theta_dot = b
    """
    psi = variational_wavefunction(x, theta)
    Hpsi = hamiltonian_action(psi, x, dx, mass, potential)
    D = tangent_vectors(x, theta, relative_step=relative_step)

    G = np.real(D.conj().T @ D) * dx
    b = np.imag(D.conj().T @ Hpsi) * dx

    velocity, _, rank, singular_values = np.linalg.lstsq(G, b, rcond=rcond)

    residual = 1j * (D @ velocity) - Hpsi
    residual_norm = np.sqrt(max(inner_product(residual, residual, dx).real, 0.0))
    zero_velocity_residual = np.sqrt(max(inner_product(Hpsi, Hpsi, dx).real, 0.0))

    return velocity, {
        "G": G,
        "b": b,
        "rank": rank,
        "singular_values": singular_values,
        "residual_norm": residual_norm,
        "zero_velocity_residual": zero_velocity_residual,
    }


def run_variational_dynamics(
    x,
    dx,
    theta0,
    mass,
    potential,
    dt=0.0005,
    steps=200,
    store_every=5,
    relative_step=1.0e-6,
    rcond=1.0e-10,
):
    """Propagate the compact real-parameter multi-Gaussian TDVP with RK4."""
    theta = np.asarray(theta0, dtype=float).copy()

    def rhs(th):
        velocity, _ = tdvp_velocity(
            x,
            dx,
            th,
            mass,
            potential,
            relative_step=relative_step,
            rcond=rcond,
        )
        return velocity

    times = []
    parameters = []
    states = []
    norms = []
    residuals = []

    def record(step):
        psi = variational_wavefunction(x, theta)
        _, info = tdvp_velocity(
            x,
            dx,
            theta,
            mass,
            potential,
            relative_step=relative_step,
            rcond=rcond,
        )

        times.append(step * dt)
        parameters.append(theta.copy())
        states.append(psi)
        norms.append((np.vdot(psi, psi) * dx).real)
        residuals.append(info["residual_norm"])

    record(0)

    for step in range(1, steps + 1):
        k1 = rhs(theta)
        k2 = rhs(theta + 0.5 * dt * k1)
        k3 = rhs(theta + 0.5 * dt * k2)
        k4 = rhs(theta + dt * k3)

        theta = theta + dt * (k1 + 2*k2 + 2*k3 + k4) / 6.0

        if step % store_every == 0:
            record(step)

    return {
        "time": np.asarray(times),
        "theta": np.asarray(parameters),
        "psi": np.asarray(states),
        "norm": np.asarray(norms),
        "residual": np.asarray(residuals),
    }
