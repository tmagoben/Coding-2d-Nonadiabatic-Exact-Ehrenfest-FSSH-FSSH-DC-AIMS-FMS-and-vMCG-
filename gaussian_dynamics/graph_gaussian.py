from dataclasses import dataclass
import numpy as np

from .gaussian_nd import analytic_overlap_equal_width
from .local_gaussian_nd import kinetic_matrix_element_equal_width


@dataclass
class GraphGaussianTBF:
    """Gaussian TBF carrying an electronic coefficient vector in its node-local frame."""
    node: object
    q: np.ndarray
    p: np.ndarray
    A: np.ndarray
    electronic_coefficients: np.ndarray

    def __post_init__(self):
        self.q = np.asarray(self.q, dtype=float)
        self.p = np.asarray(self.p, dtype=float)
        self.A = np.asarray(self.A, dtype=float)
        self.electronic_coefficients = np.asarray(
            self.electronic_coefficients, dtype=complex
        )
        if self.q.shape != self.p.shape:
            raise ValueError("q and p must have the same shape")
        if self.A.shape != (len(self.q), len(self.q)):
            raise ValueError("A has incompatible shape")


def pair_overlap_and_hamiltonian(tbf_i, tbf_j, registry, reference_node, mass_matrix):
    """Discrete local-diabatic Gaussian pair approximation.

    Electronic states at the two TBF nodes are parallel transported to one common
    reference/centroid node.  In that common gauge,

        S_ij ~= <g_i|g_j> <e_i|e_j>
        H_ij ~= <T_i|g_j> <e_i|e_j>
                + <g_i|g_j> <e_i|H_e(reference)|e_j>.

    This is a deliberately named discrete-overlap approximation.  It avoids any
    gauge-dependent direct comparison of electronic vectors living at different
    nodes, but it is not the full continuous AIMS kinetic-coupling integral.
    """
    if not np.allclose(tbf_i.A, tbf_j.A, atol=1e-12):
        raise ValueError("pair approximation currently requires equal widths")

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

    S = S_nuc * factors["overlap"]
    H = T_nuc * factors["overlap"] + S_nuc * factors["potential"]

    return S, H, factors


def build_static_graph_gaussian_matrices(
    basis,
    registry,
    mass_matrix,
    reference_selector,
):
    """Build static S,H matrices using one symmetric pair-reference choice.

    reference_selector(i,j) must return the same node for (i,j) and (j,i) if exact
    Hermiticity is desired.
    """
    n = len(basis)
    S = np.zeros((n, n), dtype=complex)
    H = np.zeros((n, n), dtype=complex)

    for i in range(n):
        for j in range(n):
            ref = reference_selector(i, j)
            S[i, j], H[i, j], _ = pair_overlap_and_hamiltonian(
                basis[i], basis[j], registry, ref, mass_matrix
            )
    return S, H


def generalized_cayley_step(coefficients, S, H, dt):
    """One norm-stable Crank-Nicolson/Cayley step for static nonorthogonal S,H."""
    C = np.asarray(coefficients, dtype=complex)
    lhs = S + 0.5j * dt * H
    rhs = (S - 0.5j * dt * H) @ C
    return np.linalg.solve(lhs, rhs)


def generalized_norm(coefficients, S):
    C = np.asarray(coefficients, dtype=complex)
    return float(np.real(np.vdot(C, S @ C)))
