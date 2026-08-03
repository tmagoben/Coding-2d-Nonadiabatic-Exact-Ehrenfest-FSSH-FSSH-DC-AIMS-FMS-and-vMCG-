"""Gaussian nuclear wavepackets and electronic-state preparation."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray


ComplexArray = NDArray[np.complex128]


def _normalize(psi: ComplexArray, volume_element: float) -> ComplexArray:
    norm_squared = volume_element * np.sum(np.abs(psi) ** 2)
    if norm_squared <= 0.0:
        raise ValueError("Cannot normalize a zero wavefunction.")
    return psi / np.sqrt(norm_squared)


def gaussian_wavepacket_1d(
    x: ArrayLike,
    *,
    center: float,
    width: float,
    momentum: float,
    dx: float,
) -> ComplexArray:
    r"""
    Return a normalized minimum-form Gaussian

        g(x) proportional to
        exp[-(x-x0)^2/(4 sigma^2) + i p0 (x-x0)].

    With this convention, |g|^2 has variance sigma^2.
    """
    x = np.asarray(x, dtype=float)
    if width <= 0.0:
        raise ValueError("width must be positive.")

    displacement = x - center
    packet = np.exp(
        -(displacement**2) / (4.0 * width**2)
        + 1j * momentum * displacement
    )
    return _normalize(packet.astype(complex), dx)


def gaussian_wavepacket_2d(
    x: ArrayLike,
    y: ArrayLike,
    *,
    center_x: float,
    center_y: float,
    width_x: float,
    width_y: float,
    momentum_x: float,
    momentum_y: float,
    dx: float,
    dy: float,
) -> ComplexArray:
    """Return a normalized separable two-dimensional Gaussian packet."""
    x, y = np.broadcast_arrays(
        np.asarray(x, dtype=float),
        np.asarray(y, dtype=float),
    )
    if width_x <= 0.0 or width_y <= 0.0:
        raise ValueError("Widths must be positive.")

    packet = np.exp(
        -(x - center_x) ** 2 / (4.0 * width_x**2)
        -(y - center_y) ** 2 / (4.0 * width_y**2)
        + 1j * momentum_x * (x - center_x)
        + 1j * momentum_y * (y - center_y)
    )
    return _normalize(packet.astype(complex), dx * dy)


def prepare_adiabatic_wavepacket_1d(
    scalar_packet: ArrayLike,
    adiabatic_vectors: ArrayLike,
    *,
    state: int,
    dx: float,
) -> ComplexArray:
    r"""
    Embed a scalar packet on a local adiabatic electronic state.

        Psi_a(x) = U_ai(x) g(x).
    """
    g = np.asarray(scalar_packet, dtype=complex)
    u = np.asarray(adiabatic_vectors, dtype=float)
    if u.shape[0] != g.size or state not in range(u.shape[-1]):
        raise ValueError("Wavepacket and adiabatic-vector dimensions differ.")

    psi = np.einsum("xa,x->ax", u[:, :, state], g, optimize=True)
    return _normalize(psi, dx)


def prepare_adiabatic_wavepacket_2d(
    scalar_packet: ArrayLike,
    adiabatic_vectors: ArrayLike,
    *,
    state: int,
    dx: float,
    dy: float,
) -> ComplexArray:
    """Embed a scalar 2D packet on one local adiabatic state."""
    g = np.asarray(scalar_packet, dtype=complex)
    u = np.asarray(adiabatic_vectors, dtype=float)
    if u.shape[:2] != g.shape or state not in range(u.shape[-1]):
        raise ValueError("Wavepacket and adiabatic-vector dimensions differ.")

    psi = np.einsum("yxa,yx->ayx", u[..., :, state], g, optimize=True)
    return _normalize(psi, dx * dy)
