"""Arbitrary-state, arbitrary-dimension initialization in a fixed electronic frame."""

from dataclasses import dataclass
import numpy as np

from .gaussian_general import gaussian_overlap_general
from .gaussian_nd import gaussian_nd
from .matrix_invariants_v213 import isometry_residual_v213, require_residual_v213


@dataclass(frozen=True)
class GridProjectionResultV213:
    coefficients: np.ndarray
    projected_wavefunction: np.ndarray
    metric: np.ndarray
    target_norm: float
    projected_norm: float
    residual_norm: float
    relative_residual: float
    fidelity: float
    condition_number: float
    nstate: int
    nuclear_dimension: int


def block_metric_fixed_frame_v213(basis, nstate):
    basis = list(basis)
    s = int(nstate)
    if not basis or s < 1:
        raise ValueError("basis and nstate must be nonempty.")
    nuclear = np.empty((len(basis), len(basis)), dtype=complex)
    for i, bi in enumerate(basis):
        for j, bj in enumerate(basis):
            nuclear[i, j] = gaussian_overlap_general(
                bi.q, bi.p, bi.A, bj.q, bj.p, bj.A
            )
    return np.kron(nuclear, np.eye(s, dtype=complex))


def transform_electronic_vector_to_local_frame_v213(vector, local_frame):
    """If ``|global> = frame |local>``, return ``c_local = frame^dagger c_global``."""
    c = np.asarray(vector, dtype=complex)
    U = np.asarray(local_frame, dtype=complex)
    if U.ndim != 2 or U.shape[0] != U.shape[1] or c.shape != (U.shape[0],):
        raise ValueError("electronic vector and local frame have incompatible shapes.")
    require_residual_v213(
        "local electronic frame unitarity", isometry_residual_v213(U), 1.0e-10
    )
    return U.conj().T @ c


def initialize_separable_block_state_v213(
    nuclear_coefficients,
    electronic_vector,
    metric,
):
    """Build and metric-normalize ``C_(i,a)=nuclear_i*electronic_a``."""
    nuclear = np.asarray(nuclear_coefficients, dtype=complex)
    electronic = np.asarray(electronic_vector, dtype=complex)
    if nuclear.ndim != 1 or electronic.ndim != 1:
        raise ValueError("nuclear and electronic coefficients must be vectors.")
    electronic_norm = float(np.real(np.vdot(electronic, electronic)))
    if electronic_norm <= 0.0:
        raise ValueError("electronic vector must be nonzero.")
    electronic = electronic / np.sqrt(electronic_norm)
    C = np.kron(nuclear, electronic)
    if getattr(metric, "shape", None) != (len(C), len(C)):
        raise ValueError("metric is incompatible with the block coefficients.")
    norm = float(np.real(np.vdot(C, metric @ C)))
    if norm <= 0.0:
        raise ValueError("separable initial state has non-positive metric norm.")
    return C / np.sqrt(norm)


def _quadrature_weights(grid_shape, weights):
    w = np.asarray(weights, dtype=float)
    if w.ndim == 0:
        w = np.full(grid_shape, float(w), dtype=float)
    else:
        try:
            w = np.broadcast_to(w, grid_shape).astype(float, copy=False)
        except ValueError as exc:
            raise ValueError("quadrature weights do not broadcast to the nuclear grid.") from exc
    if not np.all(np.isfinite(w)) or np.any(w <= 0.0):
        raise ValueError("quadrature weights must be finite and positive.")
    return w


def project_grid_wavefunction_fixed_frame_v213(
    psi_target,
    points,
    quadrature_weights,
    basis,
    *,
    rcond=1.0e-12,
):
    """Project ``psi(...,a)`` onto ``g_i(q)|a>`` for arbitrary ``d`` and ``s``.

    The electronic basis is a fixed global frame.  A later moving-frame calculation
    may transform the resulting block coefficients with the explicitly supplied frame
    matrices; no state index or spin label is assumed here.
    """
    psi = np.asarray(psi_target, dtype=complex)
    points = np.asarray(points, dtype=float)
    basis = list(basis)
    if float(rcond) <= 0.0:
        raise ValueError("rcond must be positive.")
    if points.ndim < 2:
        raise ValueError("points must have shape grid_shape+(nuclear_dimension,).")
    grid_shape = points.shape[:-1]
    d = int(points.shape[-1])
    if psi.ndim != len(grid_shape) + 1 or psi.shape[:-1] != grid_shape:
        raise ValueError("psi_target must have shape grid_shape+(nstate,).")
    s = int(psi.shape[-1])
    if s < 1 or not basis:
        raise ValueError("projection requires states and Gaussian basis functions.")
    if not np.all(np.isfinite(psi)) or not np.all(np.isfinite(points)):
        raise ValueError("projection inputs contain non-finite data.")
    for b in basis:
        if np.asarray(b.q).shape != (d,):
            raise ValueError("Gaussian nuclear dimension differs from the grid.")
    weights = _quadrature_weights(grid_shape, quadrature_weights)

    metric = block_metric_fixed_frame_v213(basis, s)
    rhs = np.zeros(len(basis) * s, dtype=complex)
    gaussians = []
    for i, b in enumerate(basis):
        g = gaussian_nd(points, b.q, b.p, b.A)
        gaussians.append(g)
        rhs[s * i : s * (i + 1)] = np.sum(
            np.conj(g)[..., None] * psi * weights[..., None],
            axis=tuple(range(len(grid_shape))),
        )

    coefficients, _, _, _ = np.linalg.lstsq(metric, rhs, rcond=float(rcond))
    projected = np.zeros_like(psi)
    for i, g in enumerate(gaussians):
        projected += g[..., None] * coefficients[s * i : s * (i + 1)]

    target_norm = float(np.real(np.sum(np.abs(psi) ** 2 * weights[..., None])))
    if target_norm <= 0.0:
        raise ValueError("psi_target must have positive quadrature norm.")
    projected_norm = float(
        np.real(np.sum(np.abs(projected) ** 2 * weights[..., None]))
    )
    residual_norm = float(
        np.real(np.sum(np.abs(projected - psi) ** 2 * weights[..., None]))
    )
    overlap = np.sum(np.conj(psi) * projected * weights[..., None])
    fidelity = float(
        abs(overlap) ** 2 / max(target_norm * projected_norm, 1.0e-30)
    )
    fidelity = min(max(fidelity, 0.0), 1.0)
    return GridProjectionResultV213(
        coefficients=coefficients,
        projected_wavefunction=projected,
        metric=metric,
        target_norm=target_norm,
        projected_norm=projected_norm,
        residual_norm=residual_norm,
        relative_residual=residual_norm / max(target_norm, 1.0e-30),
        fidelity=fidelity,
        condition_number=float(np.linalg.cond(metric)),
        nstate=s,
        nuclear_dimension=d,
    )
