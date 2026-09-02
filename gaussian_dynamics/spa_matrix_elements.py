from dataclasses import dataclass
import numpy as np

from .gaussian_nd import analytic_overlap_equal_width
from .local_gaussian_nd import (
    overlap_centroid_equal_width,
    kinetic_matrix_element_equal_width,
)


@dataclass(frozen=True)
class SPAResult:
    order: int
    overlap: complex
    kinetic: complex
    potential_zeroth: complex
    potential_first: complex
    total: complex
    saddle_point: np.ndarray
    complex_overlap_centroid: np.ndarray


def real_saddle_point_equal_width(qi, qj):
    """Maximum of |g_i g_j| for equal-width frozen Gaussians."""
    qi = np.asarray(qi, dtype=float)
    qj = np.asarray(qj, dtype=float)
    if qi.shape != qj.shape:
        raise ValueError("qi and qj must have equal shape")
    return 0.5 * (qi + qj)


def scalar_spa_matrix_element(qi, pi, qj, pj, A, value_at_saddle, gradient_at_saddle=None, order=0):
    """Gaussian matrix element of a smooth scalar quantity in SPA0/SPA1.

    For equal-width Gaussians and a Taylor expansion about the real saddle point q_c,

        f(q) ~= f(q_c) + grad f(q_c) . (q-q_c),

    the first moment is exact analytically through the complex overlap centroid mu.
    """
    if order not in (0, 1):
        raise ValueError("order must be 0 or 1")

    S = analytic_overlap_equal_width(qi, pi, qj, pj, A)
    qc = real_saddle_point_equal_width(qi, qj)
    mu = overlap_centroid_equal_width(qi, pi, qj, pj, A)

    zeroth = complex(value_at_saddle) * S
    first = 0.0 + 0.0j

    if order == 1:
        if gradient_at_saddle is None:
            raise ValueError("SPA1 requires gradient_at_saddle")
        grad = np.asarray(gradient_at_saddle, dtype=complex)
        if grad.shape != qc.shape:
            raise ValueError("gradient_at_saddle has incompatible shape")
        first = (grad @ (mu - qc)) * S

    return zeroth + first


def graph_pair_spa_result(tbf_i, tbf_j, registry, reference_node, mass_matrix, order=0):
    """Gauge-covariant graph-Gaussian SPA0/SPA1 pair matrix element.

    The electronic Hamiltonian and its first derivatives are transported/evaluated in
    one common graph reference frame.  This is a transparent electronic Taylor layer:

        SPA0: H_e(q) -> H_e(q_c)
        SPA1: H_e(q) -> H_e(q_c) + sum_a dH_e/dq_a (q_a-q_c,a)

    It is deliberately not labeled the complete production AIMS-SPA1 Hamiltonian,
    because higher derivative-coupling terms and method-specific AIMS approximations
    require additional structure.
    """
    if order not in (0, 1):
        raise ValueError("order must be 0 or 1")
    if not np.allclose(tbf_i.A, tbf_j.A, atol=1e-12):
        raise ValueError("SPA pair approximation currently requires equal widths")

    S_nuc = analytic_overlap_equal_width(
        tbf_i.q, tbf_i.p, tbf_j.q, tbf_j.p, tbf_i.A
    )
    T_nuc = kinetic_matrix_element_equal_width(
        tbf_i.q,
        tbf_i.p,
        tbf_j.q,
        tbf_j.p,
        tbf_i.A,
        mass_matrix,
    )

    factors = registry.pair_factors(
        tbf_i.node,
        tbf_i.electronic_coefficients,
        tbf_j.node,
        tbf_j.electronic_coefficients,
        reference_node,
    )

    qc = real_saddle_point_equal_width(tbf_i.q, tbf_j.q)
    mu = overlap_centroid_equal_width(
        tbf_i.q, tbf_i.p, tbf_j.q, tbf_j.p, tbf_i.A
    )

    potential0 = S_nuc * factors["potential"]
    potential1 = 0.0 + 0.0j
    if order == 1:
        deriv = np.asarray(factors["derivative_hamiltonian"], dtype=complex)
        if deriv.shape != qc.shape:
            raise ValueError("electronic derivative field does not match nuclear dimension")
        potential1 = S_nuc * (deriv @ (mu - qc))

    overlap = S_nuc * factors["overlap"]
    kinetic = T_nuc * factors["overlap"]
    total = kinetic + potential0 + potential1

    return SPAResult(
        order=order,
        overlap=overlap,
        kinetic=kinetic,
        potential_zeroth=potential0,
        potential_first=potential1,
        total=total,
        saddle_point=qc,
        complex_overlap_centroid=mu,
    )


def build_graph_gaussian_matrices_spa(basis, registry, mass_matrix, reference_selector, order=0):
    n = len(basis)
    S = np.zeros((n, n), dtype=complex)
    H = np.zeros((n, n), dtype=complex)

    for i in range(n):
        for j in range(n):
            ref = reference_selector(i, j)
            result = graph_pair_spa_result(
                basis[i], basis[j], registry, ref, mass_matrix, order=order
            )
            S[i, j] = result.overlap
            H[i, j] = result.total

    return S, H


def spa1_correction_norm(H0, H1):
    H0 = np.asarray(H0, dtype=complex)
    H1 = np.asarray(H1, dtype=complex)
    denom = max(np.linalg.norm(H1, ord="fro"), 1e-30)
    return float(np.linalg.norm(H1 - H0, ord="fro") / denom)
