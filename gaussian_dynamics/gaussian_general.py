import numpy as np


def validate_spd(A, name="A"):
    A = np.asarray(A, dtype=float)
    if A.ndim != 2 or A.shape[0] != A.shape[1]:
        raise ValueError(f"{name} must be square.")
    if not np.allclose(A, A.T, atol=1e-12):
        raise ValueError(f"{name} must be symmetric.")
    if np.min(np.linalg.eigvalsh(A)) <= 0.0:
        raise ValueError(f"{name} must be positive definite.")
    return A


def gaussian_overlap_general(qi, pi, Ai, qj, pj, Aj):
    """Exact overlap of normalized zero-chirp Gaussians with unequal widths.

    g_k(q) = N_k exp[-1/2 (q-q_k)^T A_k (q-q_k)
                     + i p_k^T(q-q_k)]

    All width matrices are real symmetric positive definite.
    """
    qi = np.asarray(qi, dtype=float)
    pi = np.asarray(pi, dtype=float)
    qj = np.asarray(qj, dtype=float)
    pj = np.asarray(pj, dtype=float)
    Ai = validate_spd(Ai, "Ai")
    Aj = validate_spd(Aj, "Aj")

    if qi.shape != qj.shape or qi.shape != pi.shape or qi.shape != pj.shape:
        raise ValueError("q and p vectors must have equal dimensions.")
    D = len(qi)
    if Ai.shape != (D, D) or Aj.shape != (D, D):
        raise ValueError("width dimensions are incompatible.")

    B = Ai + Aj
    l = Ai @ qi + Aj @ qj + 1j * (pj - pi)

    c = (
        -0.5 * qi @ Ai @ qi
        -0.5 * qj @ Aj @ qj
        +1j * pi @ qi
        -1j * pj @ qj
    )

    log_prefactor = (
        0.25 * np.linalg.slogdet(Ai)[1]
        +0.25 * np.linalg.slogdet(Aj)[1]
        -0.5 * np.linalg.slogdet(B)[1]
        +0.5 * D * np.log(2.0)
    )

    exponent = c + 0.5 * l @ np.linalg.solve(B, l)

    return np.exp(log_prefactor + exponent)


def gaussian_cross_centroid(qi, pi, Ai, qj, pj, Aj):
    """Complex centroid of the normalized cross density g_i^* g_j."""
    qi = np.asarray(qi, dtype=float)
    pi = np.asarray(pi, dtype=float)
    qj = np.asarray(qj, dtype=float)
    pj = np.asarray(pj, dtype=float)
    Ai = validate_spd(Ai, "Ai")
    Aj = validate_spd(Aj, "Aj")

    B = Ai + Aj
    l = Ai @ qi + Aj @ qj + 1j*(pj-pi)
    return np.linalg.solve(B, l)


def gaussian_cross_covariance(Ai, Aj):
    """Second central moment of the complex cross Gaussian."""
    Ai = validate_spd(Ai, "Ai")
    Aj = validate_spd(Aj, "Aj")
    return np.linalg.inv(Ai + Aj)


def real_overlap_saddle_point(qi, Ai, qj, Aj):
    """Maximum of |g_i(q) g_j(q)| for unequal-width real Gaussians."""
    qi = np.asarray(qi, dtype=float)
    qj = np.asarray(qj, dtype=float)
    Ai = validate_spd(Ai, "Ai")
    Aj = validate_spd(Aj, "Aj")
    return np.linalg.solve(Ai+Aj, Ai@qi + Aj@qj)


def gradient_matrix_element_general(qi, pi, Ai, qj, pj, Aj):
    """Exact <g_i | grad g_j> for unequal real width matrices."""
    S = gaussian_overlap_general(qi, pi, Ai, qj, pj, Aj)
    mu = gaussian_cross_centroid(qi, pi, Ai, qj, pj, Aj)
    Aj = validate_spd(Aj, "Aj")
    qj = np.asarray(qj, dtype=float)
    pj = np.asarray(pj, dtype=float)

    return (-Aj @ (mu-qj) + 1j*pj) * S


def kinetic_matrix_element_general(
    qi, pi, Ai, qj, pj, Aj, mass_matrix
):
    """Exact <g_i| -1/2 grad^T M^-1 grad |g_j> for unequal widths."""
    Ai = validate_spd(Ai, "Ai")
    Aj = validate_spd(Aj, "Aj")
    M = validate_spd(mass_matrix, "mass_matrix")
    Minv = np.linalg.inv(M)

    S = gaussian_overlap_general(qi, pi, Ai, qj, pj, Aj)
    mu = gaussian_cross_centroid(qi, pi, Ai, qj, pj, Aj)
    Sigma = gaussian_cross_covariance(Ai, Aj)

    qi = np.asarray(qi, dtype=float)
    pi = np.asarray(pi, dtype=float)
    qj = np.asarray(qj, dtype=float)
    pj = np.asarray(pj, dtype=float)

    ui = -Ai @ (mu-qi) - 1j*pi
    uj = -Aj @ (mu-qj) + 1j*pj

    fluctuation = np.trace(Ai @ Minv @ Aj @ Sigma)

    return 0.5 * S * (ui @ Minv @ uj + fluctuation)


def basis_time_matrix_element_general(
    qi,
    pi,
    Ai,
    qj,
    pj,
    Aj,
    qdot_j,
    pdot_j,
    Adot_j=None,
):
    """Exact <g_i|dot g_j> for real time-dependent zero-chirp widths.

    If Adot_j is omitted, the formula reduces to a moving frozen Gaussian.
    """
    Ai = validate_spd(Ai, "Ai")
    Aj = validate_spd(Aj, "Aj")
    qdot_j = np.asarray(qdot_j, dtype=float)
    pdot_j = np.asarray(pdot_j, dtype=float)
    qj = np.asarray(qj, dtype=float)
    pj = np.asarray(pj, dtype=float)

    if Adot_j is None:
        Adot_j = np.zeros_like(Aj)
    Adot_j = np.asarray(Adot_j, dtype=float)
    if Adot_j.shape != Aj.shape or not np.allclose(Adot_j, Adot_j.T, atol=1e-12):
        raise ValueError("Adot_j must be symmetric and match Aj.")

    S = gaussian_overlap_general(qi, pi, Ai, qj, pj, Aj)
    mu = gaussian_cross_centroid(qi, pi, Ai, qj, pj, Aj)
    Sigma = gaussian_cross_covariance(Ai, Aj)
    y = mu-qj

    normalization = 0.25*np.trace(np.linalg.solve(Aj, Adot_j))
    center = (Aj@y - 1j*pj) @ qdot_j
    momentum = 1j*y @ pdot_j
    width = -0.5*(y @ Adot_j @ y + np.trace(Adot_j @ Sigma))

    return S*(normalization + center + momentum + width)


def width_scaled(A, scale):
    A = validate_spd(A)
    scale = float(scale)
    if scale <= 0.0:
        raise ValueError("width scale must be positive.")
    return scale*A
