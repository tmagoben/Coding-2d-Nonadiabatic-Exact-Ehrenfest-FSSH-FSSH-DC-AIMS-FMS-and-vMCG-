"""Pathwise adiabatic-to-diabatic transformations."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.linalg import expm


FloatArray = NDArray[np.float64]


def integrate_adt_path(
    path_coordinate: ArrayLike,
    derivative_coupling_along_path: ArrayLike,
    initial_transformation: ArrayLike,
) -> FloatArray:
    r"""
    Integrate dA/ds = -tau_s A using midpoint matrix exponentials.
    """
    s = np.asarray(path_coordinate, dtype=float)
    tau = np.asarray(derivative_coupling_along_path, dtype=float)
    initial = np.asarray(initial_transformation, dtype=float)

    if s.ndim != 1 or tau.ndim != 3 or tau.shape[0] != s.size:
        raise ValueError("Invalid path or derivative-coupling array.")
    if np.any(np.diff(s) <= 0.0):
        raise ValueError("Path coordinates must be strictly increasing.")
    if not np.all(np.isfinite(tau)):
        raise ValueError("The path includes undefined derivative couplings.")

    transformations = np.empty_like(tau)
    transformations[0] = initial

    for p in range(s.size - 1):
        ds = s[p + 1] - s[p]
        tau_mid = 0.5 * (tau[p] + tau[p + 1])
        tau_mid = 0.5 * (tau_mid - tau_mid.T)
        candidate = expm(-tau_mid * ds) @ transformations[p]

        left, _, right_t = np.linalg.svd(candidate)
        transformations[p + 1] = left @ right_t

    return transformations


def transform_adiabatic_to_diabatic(
    energies: ArrayLike,
    transformations: ArrayLike,
) -> FloatArray:
    r"""Return V_d = A^T diag(E_ad) A at each path point."""
    e = np.asarray(energies, dtype=float)
    a = np.asarray(transformations, dtype=float)
    if e.ndim != 2 or a.shape != (e.shape[0], e.shape[1], e.shape[1]):
        raise ValueError("Incompatible energies and transformations.")

    result = np.empty_like(a)
    for p in range(e.shape[0]):
        result[p] = a[p].T @ np.diag(e[p]) @ a[p]
    return result
