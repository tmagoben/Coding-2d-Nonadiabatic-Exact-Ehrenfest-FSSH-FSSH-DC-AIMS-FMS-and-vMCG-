"""State tracking and nonadiabatic derivative couplings."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.optimize import linear_sum_assignment


FloatArray = NDArray[np.float64]


def diagonalize_path(
    hamiltonians: ArrayLike,
) -> tuple[FloatArray, FloatArray]:
    r"""
    Diagonalize a sequence of real symmetric Hamiltonians and track states.

    State assignment maximizes the absolute overlap with the previous point.
    Each state is then sign-corrected to have positive overlap with its
    predecessor.
    """
    h = np.asarray(hamiltonians, dtype=float)
    if h.ndim != 3 or h.shape[1] != h.shape[2]:
        raise ValueError("Expected shape (n_points,n_states,n_states).")
    if not np.allclose(h, np.swapaxes(h, -1, -2), atol=1.0e-12):
        raise ValueError("Hamiltonians must be real symmetric.")

    n_points, n_states, _ = h.shape
    energies = np.empty((n_points, n_states), dtype=float)
    vectors = np.empty((n_points, n_states, n_states), dtype=float)

    energies[0], vectors[0] = np.linalg.eigh(h[0])

    for p in range(1, n_points):
        e_now, u_now = np.linalg.eigh(h[p])
        overlap = vectors[p - 1].T @ u_now
        rows, columns = linear_sum_assignment(-np.abs(overlap))

        permutation = np.empty(n_states, dtype=int)
        permutation[rows] = columns
        e_now = e_now[permutation]
        u_now = u_now[:, permutation]

        signs = np.where(
            np.diag(vectors[p - 1].T @ u_now) < 0.0,
            -1.0,
            1.0,
        )
        vectors[p] = u_now * signs[np.newaxis, :]
        energies[p] = e_now

    return energies, vectors


def hellmann_feynman_derivative_couplings(
    energies: ArrayLike,
    eigenvectors: ArrayLike,
    hamiltonian_derivatives: ArrayLike,
    *,
    gap_tolerance: float = 1.0e-12,
) -> FloatArray:
    r"""
    Compute off-diagonal derivative couplings

        tau_ij^alpha =
        <phi_i|dH/dR_alpha|phi_j> / (E_j-E_i).

    The diagonal gauge terms are set to zero and the real coupling matrix is
    explicitly antisymmetrized.
    """
    e = np.asarray(energies, dtype=float)
    u = np.asarray(eigenvectors, dtype=float)
    dh = np.asarray(hamiltonian_derivatives, dtype=float)

    if e.ndim != 2 or u.ndim != 3 or dh.ndim != 4:
        raise ValueError("Incompatible energy, vector, or derivative arrays.")

    n_points, n_states = e.shape
    n_coordinates = dh.shape[1]
    tau = np.zeros(
        (n_points, n_coordinates, n_states, n_states),
        dtype=float,
    )

    for p in range(n_points):
        for alpha in range(n_coordinates):
            derivative_ad = u[p].T @ dh[p, alpha] @ u[p]
            for i in range(n_states):
                for j in range(i + 1, n_states):
                    gap = e[p, j] - e[p, i]
                    value = (
                        np.nan
                        if abs(gap) <= gap_tolerance
                        else derivative_ad[i, j] / gap
                    )
                    tau[p, alpha, i, j] = value
                    tau[p, alpha, j, i] = -value
    return tau


def finite_difference_derivative_couplings_1d(
    coordinate: ArrayLike,
    eigenvectors: ArrayLike,
) -> FloatArray:
    r"""Compute tau = U^T dU/dx from phase-tracked eigenvectors."""
    x = np.asarray(coordinate, dtype=float)
    u = np.asarray(eigenvectors, dtype=float)
    if x.ndim != 1 or u.ndim != 3 or u.shape[0] != x.size:
        raise ValueError("Incompatible coordinate and eigenvector arrays.")

    du_dx = np.gradient(u, x, axis=0, edge_order=2)
    tau = np.einsum("pbi,pbj->pij", u, du_dx, optimize=True)
    return 0.5 * (tau - np.swapaxes(tau, -1, -2))
