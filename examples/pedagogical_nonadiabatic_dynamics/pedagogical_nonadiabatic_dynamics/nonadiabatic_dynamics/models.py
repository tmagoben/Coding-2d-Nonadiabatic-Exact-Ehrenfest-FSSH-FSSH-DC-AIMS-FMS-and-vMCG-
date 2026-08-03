"""Analytic two-state model Hamiltonians."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray


FloatArray = NDArray[np.float64]


def smooth_single_avoided_crossing(
    x: ArrayLike,
    *,
    amplitude: float = 0.010,
    slope_scale: float = 1.600,
    coupling: float = 0.005,
    coupling_decay: float = 1.000,
) -> tuple[FloatArray, FloatArray]:
    r"""
    Return a fully smooth Tully-like single avoided-crossing model.

    The diabatic matrix is

        V_11(x) = A tanh(Bx),
        V_22(x) = -V_11(x),
        V_12(x) = C exp(-D x^2).

    This retains the familiar scattering topology of the Tully single
    avoided crossing while making every derivative smooth at x=0, which is
    useful for a finite-difference derivative-coupling tutorial.
    """
    x = np.asarray(x, dtype=float)
    if x.ndim != 1:
        raise ValueError("x must be one-dimensional.")

    v11 = amplitude * np.tanh(slope_scale * x)
    v12 = coupling * np.exp(-coupling_decay * x**2)

    potential = np.zeros((x.size, 2, 2), dtype=float)
    potential[:, 0, 0] = v11
    potential[:, 1, 1] = -v11
    potential[:, 0, 1] = v12
    potential[:, 1, 0] = v12

    sech_squared = 1.0 / np.cosh(slope_scale * x) ** 2
    dv11 = amplitude * slope_scale * sech_squared
    dv12 = -2.0 * coupling_decay * x * v12

    derivative = np.zeros((x.size, 1, 2, 2), dtype=float)
    derivative[:, 0, 0, 0] = dv11
    derivative[:, 0, 1, 1] = -dv11
    derivative[:, 0, 0, 1] = dv12
    derivative[:, 0, 1, 0] = dv12

    return potential, derivative


# Backward-compatible name for early versions of this tutorial.
tully_single_avoided_crossing = smooth_single_avoided_crossing


def linear_vibronic_coupling_2d(
    x: ArrayLike,
    y: ArrayLike,
    *,
    kappa: float = 0.020,
    lambda_: float = 0.015,
    common_curvature: float = 0.002,
) -> tuple[FloatArray, FloatArray]:
    r"""
    Return a two-state linear vibronic-coupling conical-intersection model.

        V_d(x,y) = 1/2 k (x^2+y^2) I
                   + [[kappa*x, lambda*y],
                      [lambda*y, -kappa*x]].

    Returns
    -------
    potential
        Shape ``broadcast(x,y).shape + (2,2)``.
    derivatives
        Shape ``broadcast(x,y).shape + (2,2,2)`` where coordinate index 0
        is dV/dx and index 1 is dV/dy.
    """
    x, y = np.broadcast_arrays(
        np.asarray(x, dtype=float),
        np.asarray(y, dtype=float),
    )
    shape = x.shape
    common = 0.5 * common_curvature * (x**2 + y**2)

    potential = np.zeros(shape + (2, 2), dtype=float)
    potential[..., 0, 0] = common + kappa * x
    potential[..., 1, 1] = common - kappa * x
    potential[..., 0, 1] = lambda_ * y
    potential[..., 1, 0] = lambda_ * y

    derivatives = np.zeros(shape + (2, 2, 2), dtype=float)

    derivatives[..., 0, 0, 0] = common_curvature * x + kappa
    derivatives[..., 0, 1, 1] = common_curvature * x - kappa

    derivatives[..., 1, 0, 0] = common_curvature * y
    derivatives[..., 1, 1, 1] = common_curvature * y
    derivatives[..., 1, 0, 1] = lambda_
    derivatives[..., 1, 1, 0] = lambda_

    return potential, derivatives


def lvc_analytic_eigensystem(
    x: ArrayLike,
    y: ArrayLike,
    *,
    kappa: float = 0.020,
    lambda_: float = 0.015,
    common_curvature: float = 0.002,
) -> tuple[FloatArray, FloatArray]:
    """Return a deterministic local real adiabatic gauge for the LVC model."""
    x, y = np.broadcast_arrays(
        np.asarray(x, dtype=float),
        np.asarray(y, dtype=float),
    )
    common = 0.5 * common_curvature * (x**2 + y**2)
    qx = kappa * x
    qy = lambda_ * y
    radius = np.sqrt(qx**2 + qy**2)

    energies = np.stack((common - radius, common + radius), axis=-1)

    theta = np.arctan2(qy, qx)
    beta = 0.5 * theta
    sin_beta = np.sin(beta)
    cos_beta = np.cos(beta)

    vectors = np.zeros(x.shape + (2, 2), dtype=float)
    vectors[..., 0, 0] = -sin_beta
    vectors[..., 1, 0] = cos_beta
    vectors[..., 0, 1] = cos_beta
    vectors[..., 1, 1] = sin_beta
    return energies, vectors


def lvc_analytic_derivative_coupling(
    x: ArrayLike,
    y: ArrayLike,
    *,
    kappa: float = 0.020,
    lambda_: float = 0.015,
    singular_radius: float = 1.0e-12,
) -> tuple[FloatArray, FloatArray, NDArray[np.bool_]]:
    r"""
    Return tau_01 = <phi_0|grad phi_1> for the chosen local LVC gauge.

        tau_x = -kappa*lambda*y / [2(kappa^2 x^2 + lambda^2 y^2)]
        tau_y =  kappa*lambda*x / [2(kappa^2 x^2 + lambda^2 y^2)]
    """
    x, y = np.broadcast_arrays(
        np.asarray(x, dtype=float),
        np.asarray(y, dtype=float),
    )
    denominator = (kappa * x) ** 2 + (lambda_ * y) ** 2
    valid = denominator > singular_radius**2

    tau_x = np.full(x.shape, np.nan, dtype=float)
    tau_y = np.full(x.shape, np.nan, dtype=float)
    tau_x[valid] = (
        -0.5 * kappa * lambda_ * y[valid] / denominator[valid]
    )
    tau_y[valid] = (
        0.5 * kappa * lambda_ * x[valid] / denominator[valid]
    )
    return tau_x, tau_y, valid
