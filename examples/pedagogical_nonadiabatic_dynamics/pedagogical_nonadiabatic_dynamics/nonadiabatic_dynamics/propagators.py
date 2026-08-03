"""FFT split-operator propagation for coupled two-state wavefunctions."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray


ComplexArray = NDArray[np.complex128]


def potential_propagator_2x2(
    potential: ArrayLike,
    time_step: float,
) -> ComplexArray:
    r"""
    Compute exp(-i V dt) for a real symmetric 2x2 matrix at every grid point.

    Write

        V = v0 I + hx sigma_x + hz sigma_z.

    Then

        exp(-iVdt) = exp(-iv0dt)
        [cos(qdt) I - i sin(qdt)(hx sigma_x+hz sigma_z)/q],

    where q = sqrt(hx^2+hz^2).
    """
    v = np.asarray(potential, dtype=float)
    if v.shape[-2:] != (2, 2):
        raise ValueError("Potential must end in a 2x2 matrix.")
    if not np.allclose(v, np.swapaxes(v, -1, -2), atol=1.0e-12):
        raise ValueError("Potential must be real symmetric.")

    v0 = 0.5 * (v[..., 0, 0] + v[..., 1, 1])
    hz = 0.5 * (v[..., 0, 0] - v[..., 1, 1])
    hx = v[..., 0, 1]
    q = np.sqrt(hx**2 + hz**2)

    phase = np.exp(-1j * v0 * time_step)
    cosine = np.cos(q * time_step)
    sinc_factor = np.empty_like(q)
    small = q < 1.0e-14
    sinc_factor[small] = time_step
    sinc_factor[~small] = np.sin(q[~small] * time_step) / q[~small]

    propagator = np.empty(v.shape, dtype=complex)
    propagator[..., 0, 0] = phase * (
        cosine - 1j * sinc_factor * hz
    )
    propagator[..., 1, 1] = phase * (
        cosine + 1j * sinc_factor * hz
    )
    propagator[..., 0, 1] = phase * (
        -1j * sinc_factor * hx
    )
    propagator[..., 1, 0] = propagator[..., 0, 1]
    return propagator


def _apply_local_1d(
    propagator: ComplexArray,
    psi: ComplexArray,
) -> ComplexArray:
    return np.einsum("xab,bx->ax", propagator, psi, optimize=True)


def _apply_local_2d(
    propagator: ComplexArray,
    psi: ComplexArray,
) -> ComplexArray:
    return np.einsum("yxab,byx->ayx", propagator, psi, optimize=True)


def split_operator_step_1d(
    psi: ArrayLike,
    potential_half_step: ArrayLike,
    kinetic_full_step: ArrayLike,
) -> ComplexArray:
    r"""
    One second-order Strang step

        exp(-iVdt/2) exp(-iTdt) exp(-iVdt/2).
    """
    state = np.asarray(psi, dtype=complex)
    u_v = np.asarray(potential_half_step, dtype=complex)
    u_t = np.asarray(kinetic_full_step, dtype=complex)

    state = _apply_local_1d(u_v, state)
    state_k = np.fft.fft(state, axis=-1)
    state_k *= u_t[np.newaxis, :]
    state = np.fft.ifft(state_k, axis=-1)
    return _apply_local_1d(u_v, state)


def split_operator_step_2d(
    psi: ArrayLike,
    potential_half_step: ArrayLike,
    kinetic_full_step: ArrayLike,
) -> ComplexArray:
    """One two-dimensional second-order Strang split-operator step."""
    state = np.asarray(psi, dtype=complex)
    u_v = np.asarray(potential_half_step, dtype=complex)
    u_t = np.asarray(kinetic_full_step, dtype=complex)

    state = _apply_local_2d(u_v, state)
    state_k = np.fft.fftn(state, axes=(-2, -1))
    state_k *= u_t[np.newaxis, :, :]
    state = np.fft.ifftn(state_k, axes=(-2, -1))
    return _apply_local_2d(u_v, state)
