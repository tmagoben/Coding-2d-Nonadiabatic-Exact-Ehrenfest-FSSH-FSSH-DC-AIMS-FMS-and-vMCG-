"""Norms, populations, and comparisons between propagated wavefunctions."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray


FloatArray = NDArray[np.float64]


def norm(psi: ArrayLike, volume_element: float) -> float:
    state = np.asarray(psi, dtype=complex)
    return float(volume_element * np.sum(np.abs(state) ** 2))


def diabatic_populations(
    psi: ArrayLike,
    volume_element: float,
) -> FloatArray:
    state = np.asarray(psi, dtype=complex)
    spatial_axes = tuple(range(1, state.ndim))
    return (
        volume_element
        * np.sum(np.abs(state) ** 2, axis=spatial_axes)
    ).real


def adiabatic_populations_1d(
    psi_diabatic: ArrayLike,
    adiabatic_vectors: ArrayLike,
    dx: float,
) -> FloatArray:
    r"""Transform Psi_ad = U^T Psi_d and integrate each state population."""
    psi = np.asarray(psi_diabatic, dtype=complex)
    u = np.asarray(adiabatic_vectors, dtype=float)
    psi_ad = np.einsum("xai,ax->ix", u, psi, optimize=True)
    return (dx * np.sum(np.abs(psi_ad) ** 2, axis=1)).real


def adiabatic_populations_2d(
    psi_diabatic: ArrayLike,
    adiabatic_vectors: ArrayLike,
    dx: float,
    dy: float,
) -> FloatArray:
    """Two-dimensional local transformation from diabatic to adiabatic."""
    psi = np.asarray(psi_diabatic, dtype=complex)
    u = np.asarray(adiabatic_vectors, dtype=float)
    psi_ad = np.einsum("yxai,ayx->iyx", u, psi, optimize=True)
    return (
        dx * dy * np.sum(np.abs(psi_ad) ** 2, axis=(1, 2))
    ).real


def phase_aligned_error(
    reference: ArrayLike,
    candidate: ArrayLike,
    volume_element: float,
) -> tuple[float, float]:
    """Return fidelity and phase-aligned L2 error."""
    ref = np.asarray(reference, dtype=complex)
    test = np.asarray(candidate, dtype=complex)
    overlap = volume_element * np.vdot(ref.reshape(-1), test.reshape(-1))
    norm_ref = volume_element * np.vdot(ref.reshape(-1), ref.reshape(-1)).real
    norm_test = volume_element * np.vdot(test.reshape(-1), test.reshape(-1)).real
    fidelity = float(abs(overlap) ** 2 / (norm_ref * norm_test))

    phase = 1.0 + 0.0j if abs(overlap) == 0.0 else overlap / abs(overlap)
    difference = test - phase * ref
    error = float(np.sqrt(volume_element * np.sum(np.abs(difference) ** 2)))
    return fidelity, error
