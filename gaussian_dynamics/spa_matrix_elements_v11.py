from dataclasses import dataclass
import numpy as np

from .gaussian_general import (
    gaussian_overlap_general,
    gaussian_cross_centroid,
    real_overlap_saddle_point,
    kinetic_matrix_element_general,
)


@dataclass(frozen=True)
class GeneralSPAResult:
    order: int
    overlap: complex
    kinetic: complex
    potential_zeroth: complex
    potential_first: complex
    total: complex
    saddle_point: np.ndarray
    complex_overlap_centroid: np.ndarray


def graph_pair_spa_result_general(
    tbf_i,
    tbf_j,
    registry,
    reference_node,
    mass_matrix,
    order=0,
):
    """Gauge-covariant SPA0/SPA1 pair element for unequal Gaussian widths."""
    if order not in (0,1):
        raise ValueError("order must be 0 or 1.")

    S_nuc = gaussian_overlap_general(
        tbf_i.q, tbf_i.p, tbf_i.A,
        tbf_j.q, tbf_j.p, tbf_j.A,
    )
    T_nuc = kinetic_matrix_element_general(
        tbf_i.q, tbf_i.p, tbf_i.A,
        tbf_j.q, tbf_j.p, tbf_j.A,
        mass_matrix,
    )

    factors = registry.pair_factors(
        tbf_i.node,
        tbf_i.electronic_coefficients,
        tbf_j.node,
        tbf_j.electronic_coefficients,
        reference_node,
    )

    qc = real_overlap_saddle_point(
        tbf_i.q, tbf_i.A, tbf_j.q, tbf_j.A
    )
    mu = gaussian_cross_centroid(
        tbf_i.q, tbf_i.p, tbf_i.A,
        tbf_j.q, tbf_j.p, tbf_j.A,
    )

    V0 = S_nuc*factors["potential"]
    V1 = 0.0+0.0j

    if order == 1:
        derivative = np.asarray(
            factors["derivative_hamiltonian"],
            dtype=complex,
        )
        if derivative.shape != qc.shape:
            raise ValueError("electronic derivative field dimension mismatch.")
        V1 = S_nuc*(derivative @ (mu-qc))

    overlap = S_nuc*factors["overlap"]
    kinetic = T_nuc*factors["overlap"]

    return GeneralSPAResult(
        order=order,
        overlap=overlap,
        kinetic=kinetic,
        potential_zeroth=V0,
        potential_first=V1,
        total=kinetic+V0+V1,
        saddle_point=qc,
        complex_overlap_centroid=mu,
    )


def build_graph_gaussian_matrices_spa_general(
    basis,
    registry,
    mass_matrix,
    reference_selector,
    order=0,
):
    n = len(basis)
    S = np.zeros((n,n),dtype=complex)
    H = np.zeros((n,n),dtype=complex)

    for i in range(n):
        for j in range(n):
            result = graph_pair_spa_result_general(
                basis[i],
                basis[j],
                registry,
                reference_selector(i,j),
                mass_matrix,
                order=order,
            )
            S[i,j] = result.overlap
            H[i,j] = result.total

    return S,H


def spa1_correction_norm_general(H0,H1):
    H0=np.asarray(H0,complex)
    H1=np.asarray(H1,complex)
    denom=max(np.linalg.norm(H1,ord="fro"),1e-30)
    return float(np.linalg.norm(H1-H0,ord="fro")/denom)
