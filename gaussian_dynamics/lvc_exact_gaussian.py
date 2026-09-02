from dataclasses import dataclass
import numpy as np

from .ci2d import LVC2DParameters
from .gaussian_general import (
    gaussian_overlap_general,
    gaussian_cross_centroid,
    gaussian_cross_covariance,
    kinetic_matrix_element_general,
    basis_time_matrix_element_general,
)


IDENTITY_2 = np.eye(2, dtype=complex)
SIGMA_X = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=complex)
SIGMA_Z = np.array([[1.0, 0.0], [0.0, -1.0]], dtype=complex)


@dataclass(frozen=True)
class ExactLVCPairResult:
    overlap: complex
    kinetic: complex
    potential: complex
    total: complex
    nuclear_overlap: complex
    electronic_overlap: complex
    potential_matrix_diabatic: np.ndarray
    complex_centroid: np.ndarray
    cross_covariance: np.ndarray


def tbf_electronic_spinor(tbf, provider):
    """Electronic spinor used by a basis function in the global diabatic basis.

    v0.12 LocalDiabaticTBF objects carry an explicit `spinor`.  Older TBF objects do
    not, so the fallback is the instantaneous adiabatic eigenvector at the TBF center.
    """
    if hasattr(tbf, "spinor"):
        v = np.asarray(tbf.spinor, dtype=complex)
        n = np.linalg.norm(v)
        if n <= 0.0:
            raise ValueError("TBF electronic spinor cannot be zero.")
        return v/n

    point = provider.evaluate(np.asarray(tbf.q, dtype=float))
    return np.asarray(point.frame[:, int(tbf.state)], dtype=complex)


def center_adiabatic_spinor(tbf, provider):
    """Backward-compatible alias for the instantaneous center adiabatic spinor."""
    point = provider.evaluate(np.asarray(tbf.q, dtype=float))
    return np.asarray(point.frame[:, int(tbf.state)], dtype=complex)


def center_spinor_time_derivative(tbf, provider, qdot=None):
    r"""Time derivative of the center adiabatic spinor in the provider's global basis.

    The provider contract is

        d_ij,alpha = <phi_i | partial_alpha phi_j>.

    Completeness gives

        partial_alpha |phi_j>
          = sum_i |phi_i> d_ij,alpha,

    so

        d/dt |phi_j>
          = Phi(R_j) [v . d]_:,j.
    """
    point = provider.evaluate(np.asarray(tbf.q, dtype=float))

    if qdot is None:
        qdot = np.linalg.solve(
            np.asarray(point.mass_matrix, dtype=float),
            np.asarray(tbf.p, dtype=float),
        )
    qdot = np.asarray(qdot, dtype=float)

    directional = np.einsum(
        "ija,a->ij",
        np.asarray(point.nac, dtype=float),
        qdot,
    )
    return np.asarray(point.frame, dtype=complex) @ directional[:, int(tbf.state)]


def exact_lvc_potential_matrix_element(
    qi,
    pi,
    Ai,
    qj,
    pj,
    Aj,
    params=LVC2DParameters(),
):
    r"""Return the 2x2 diabatic matrix <g_i | V_d(R) | g_j> exactly.

    For the repository LVC model

        V_d(R)
        = 1/2 omega^2 (x^2+y^2) I
          + kappa x sigma_z
          + lambda y sigma_x.

    Since this polynomial is at most quadratic, only first and second Gaussian cross
    moments are required.
    """
    qi = np.asarray(qi, dtype=float)
    qj = np.asarray(qj, dtype=float)

    if qi.shape != (2,) or qj.shape != (2,):
        raise ValueError("The exact LVC matrix element is specialized to 2D.")

    S = gaussian_overlap_general(qi, pi, Ai, qj, pj, Aj)
    mu = gaussian_cross_centroid(qi, pi, Ai, qj, pj, Aj)
    Sigma = gaussian_cross_covariance(Ai, Aj)

    # Cross moments are algebraic moments, not |mu|^2.
    second_sum = (
        mu[0] * mu[0]
        + mu[1] * mu[1]
        + Sigma[0, 0]
        + Sigma[1, 1]
    )

    common = 0.5 * params.omega**2 * second_sum

    matrix = S * (
        common * IDENTITY_2
        + params.kappa * mu[0] * SIGMA_Z
        + params.lam * mu[1] * SIGMA_X
    )
    return matrix


def exact_lvc_pair_result(
    tbf_i,
    tbf_j,
    provider,
):
    """Exact analytic LVC matrix element for center-frozen electronic spinors.

    Basis function:
        |Xi_i(R)> = g_i(R) |u_i>,
    where |u_i> is the adiabatic eigenvector at the TBF center represented in the
    global diabatic basis and is held constant with respect to the integration
    coordinate R for that basis function.

    This is a local-diabatic / frozen-spinor Gaussian ansatz, not the full
    coordinate-dependent Born-Huang adiabatic TBF ansatz.
    """
    point_i = provider.evaluate(np.asarray(tbf_i.q, dtype=float))
    point_j = provider.evaluate(np.asarray(tbf_j.q, dtype=float))

    if not np.allclose(point_i.mass_matrix, point_j.mass_matrix, atol=1e-12):
        raise ValueError("Exact LVC pair helper currently assumes one mass matrix.")

    ui = tbf_electronic_spinor(tbf_i, provider)
    uj = tbf_electronic_spinor(tbf_j, provider)

    S_nuc = gaussian_overlap_general(
        tbf_i.q, tbf_i.p, tbf_i.A,
        tbf_j.q, tbf_j.p, tbf_j.A,
    )
    electronic_overlap = np.vdot(ui, uj)
    overlap = S_nuc * electronic_overlap

    T_nuc = kinetic_matrix_element_general(
        tbf_i.q, tbf_i.p, tbf_i.A,
        tbf_j.q, tbf_j.p, tbf_j.A,
        point_i.mass_matrix,
    )
    kinetic = T_nuc * electronic_overlap

    params = getattr(provider, "params", LVC2DParameters())
    V_matrix = exact_lvc_potential_matrix_element(
        tbf_i.q, tbf_i.p, tbf_i.A,
        tbf_j.q, tbf_j.p, tbf_j.A,
        params=params,
    )
    potential = np.vdot(ui, V_matrix @ uj)

    return ExactLVCPairResult(
        overlap=overlap,
        kinetic=kinetic,
        potential=potential,
        total=kinetic + potential,
        nuclear_overlap=S_nuc,
        electronic_overlap=electronic_overlap,
        potential_matrix_diabatic=V_matrix,
        complex_centroid=gaussian_cross_centroid(
            tbf_i.q, tbf_i.p, tbf_i.A,
            tbf_j.q, tbf_j.p, tbf_j.A,
        ),
        cross_covariance=gaussian_cross_covariance(tbf_i.A, tbf_j.A),
    )


def build_exact_lvc_gaussian_matrices(basis, provider):
    """Build exact S and H matrices for the analytic 2D LVC frozen-spinor basis."""
    n = len(basis)
    S = np.zeros((n, n), dtype=complex)
    H = np.zeros((n, n), dtype=complex)

    for i in range(n):
        for j in range(n):
            result = exact_lvc_pair_result(basis[i], basis[j], provider)
            S[i, j] = result.overlap
            H[i, j] = result.total

    return S, H


def exact_lvc_basis_time_matrix(
    basis,
    provider,
    qdots,
    pdots,
    Adots=None,
):
    r"""Exact moving-basis matrix <Xi_i | dot Xi_j> for center-frozen spinors.

    With
        Xi_j = g_j u_j,

        <Xi_i|dot Xi_j>
          = <g_i|dot g_j> <u_i|u_j>
            + <g_i|g_j> <u_i|dot u_j>.

    The second term is the electronic contribution missing from a purely nuclear
    moving-basis seed.
    """
    n = len(basis)
    qdots = np.asarray(qdots, dtype=float)
    pdots = np.asarray(pdots, dtype=float)

    if n == 0:
        raise ValueError("basis cannot be empty.")
    nq = len(basis[0].q)

    if qdots.shape != (n, nq) or pdots.shape != (n, nq):
        raise ValueError("qdots/pdots have incompatible shape.")

    if Adots is None:
        Adots = [np.zeros_like(b.A) for b in basis]
    if len(Adots) != n:
        raise ValueError("Adots must contain one matrix per basis function.")

    spinors = [tbf_electronic_spinor(b, provider) for b in basis]
    spinor_dots = []
    for i, b in enumerate(basis):
        if hasattr(b, "spinor"):
            # Explicit stored local-diabatic spinors are taken as parallel transported
            # unless the caller supplies a different derivative separately.
            spinor_dots.append(np.zeros_like(spinors[i]))
        else:
            spinor_dots.append(
                center_spinor_time_derivative(b, provider, qdot=qdots[i])
            )

    T = np.zeros((n, n), dtype=complex)

    for i in range(n):
        for j in range(n):
            S_nuc = gaussian_overlap_general(
                basis[i].q, basis[i].p, basis[i].A,
                basis[j].q, basis[j].p, basis[j].A,
            )
            T_nuc = basis_time_matrix_element_general(
                basis[i].q,
                basis[i].p,
                basis[i].A,
                basis[j].q,
                basis[j].p,
                basis[j].A,
                qdots[j],
                pdots[j],
                Adots[j],
            )

            electronic_overlap = np.vdot(spinors[i], spinors[j])
            electronic_time = np.vdot(spinors[i], spinor_dots[j])

            T[i, j] = (
                T_nuc * electronic_overlap
                + S_nuc * electronic_time
            )

    return T


def exact_lvc_basis_time_matrix_with_spinor_derivatives(
    basis,
    provider,
    qdots,
    pdots,
    spinor_dots,
    Adots=None,
):
    """General <Xi_i|dot Xi_j> with caller-supplied electronic spinor derivatives."""
    n=len(basis)
    qdots=np.asarray(qdots,dtype=float)
    pdots=np.asarray(pdots,dtype=float)
    spinors=[tbf_electronic_spinor(b,provider) for b in basis]
    spinor_dots=[np.asarray(v,dtype=complex) for v in spinor_dots]

    if len(spinor_dots)!=n:
        raise ValueError("spinor_dots must contain one vector per TBF.")
    if Adots is None:
        Adots=[np.zeros_like(b.A) for b in basis]

    T=np.zeros((n,n),dtype=complex)

    for i in range(n):
        for j in range(n):
            S_nuc=gaussian_overlap_general(
                basis[i].q,basis[i].p,basis[i].A,
                basis[j].q,basis[j].p,basis[j].A,
            )
            T_nuc=basis_time_matrix_element_general(
                basis[i].q,basis[i].p,basis[i].A,
                basis[j].q,basis[j].p,basis[j].A,
                qdots[j],pdots[j],Adots[j],
            )
            T[i,j]=(
                T_nuc*np.vdot(spinors[i],spinors[j])
                +S_nuc*np.vdot(spinors[i],spinor_dots[j])
            )

    return T
