import numpy as np


def _validate_width(A):
    A = np.asarray(A, dtype=float)
    if A.ndim != 2 or A.shape[0] != A.shape[1]:
        raise ValueError("A must be a square matrix.")
    if not np.allclose(A, A.T, atol=1e-12):
        raise ValueError("A must be symmetric.")
    eig = np.linalg.eigvalsh(A)
    if np.min(eig) <= 0.0:
        raise ValueError("A must be positive definite.")
    return A


def gaussian_nd(points, q, p, A, K=None):
    """Normalized D-dimensional Gaussian evaluated at points[...,D]."""
    points = np.asarray(points, dtype=float)
    q = np.asarray(q, dtype=float)
    p = np.asarray(p, dtype=float)
    A = _validate_width(A)

    D = len(q)
    if points.shape[-1] != D or p.shape != (D,) or A.shape != (D, D):
        raise ValueError("Incompatible multidimensional Gaussian shapes.")

    if K is None:
        K = np.zeros_like(A)
    K = np.asarray(K, dtype=float)
    if K.shape != A.shape or not np.allclose(K, K.T, atol=1e-12):
        raise ValueError("K must be symmetric and match A.")

    xi = points - q
    real_quad = np.einsum("...i,ij,...j->...", xi, A, xi)
    phase_lin = np.einsum("...i,i->...", xi, p)
    phase_quad = np.einsum("...i,ij,...j->...", xi, K, xi)

    norm = (np.linalg.det(A) / np.pi**D) ** 0.25

    return norm * np.exp(
        -0.5*real_quad
        + 1j*phase_lin
        + 0.5j*phase_quad
    )


def gaussian_nd_gradient(points, q, p, A, K=None):
    points = np.asarray(points, dtype=float)
    A = _validate_width(A)
    if K is None:
        K = np.zeros_like(A)
    K = np.asarray(K, dtype=float)

    g = gaussian_nd(points, q, p, A, K)
    xi = points - np.asarray(q, dtype=float)
    Z = A - 1j*K
    f = -np.einsum("ij,...j->...i", Z, xi) + 1j*np.asarray(p)
    return f * g[..., None]


def gaussian_nd_laplacian(points, q, p, A, K=None):
    points = np.asarray(points, dtype=float)
    A = _validate_width(A)
    if K is None:
        K = np.zeros_like(A)
    K = np.asarray(K, dtype=float)

    g = gaussian_nd(points, q, p, A, K)
    xi = points - np.asarray(q, dtype=float)
    Z = A - 1j*K
    f = -np.einsum("ij,...j->...i", Z, xi) + 1j*np.asarray(p)
    factor = np.einsum("...i,...i->...", f, f) - np.trace(Z)
    return factor * g


def kinetic_on_gaussian_nd(points, q, p, A, mass=1.0, K=None):
    return -gaussian_nd_laplacian(points, q, p, A, K)/(2.0*mass)


def analytic_overlap_equal_width(qi, pi, qj, pj, A):
    A = _validate_width(A)
    qi = np.asarray(qi, float)
    qj = np.asarray(qj, float)
    pi = np.asarray(pi, float)
    pj = np.asarray(pj, float)

    dq = qi-qj
    dp = pi-pj
    Ainv = np.linalg.inv(A)

    return np.exp(
        -0.25*dq @ A @ dq
        -0.25*dp @ Ainv @ dp
        +0.5j*(pi+pj) @ dq
    )


def gaussian_nd_time_derivative(points, q, p, A, qdot, pdot):
    """Fixed-width, zero-chirp basis derivative from center motion."""
    g = gaussian_nd(points, q, p, A)
    xi = np.asarray(points)-np.asarray(q)
    dq_factor = np.einsum("ij,...j->...i", A, xi) - 1j*np.asarray(p)
    dp_factor = 1j*xi
    factor = (
        np.einsum("...i,i->...", dq_factor, np.asarray(qdot))
        + np.einsum("...i,i->...", dp_factor, np.asarray(pdot))
    )
    return factor*g
