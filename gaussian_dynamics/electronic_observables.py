import numpy as np

from .gaussian_general import gaussian_overlap_general


def _one_hot(state, dimension):
    c = np.zeros(dimension, dtype=complex)
    c[int(state)] = 1.0
    return c




def reduced_electronic_density_from_vectors(
    coefficients,
    basis,
    electronic_vectors,
    normalize=True,
):
    """Reduced electronic density from nuclear Gaussians and common-basis vectors."""
    C = np.asarray(coefficients, dtype=complex)
    vectors = [np.asarray(v, dtype=complex) for v in electronic_vectors]

    if len(C) != len(basis) or len(vectors) != len(basis):
        raise ValueError("coefficients, basis, and electronic_vectors must match.")

    dim = len(vectors[0])
    if any(v.shape != (dim,) for v in vectors):
        raise ValueError("electronic vectors must have one common dimension.")

    rho = np.zeros((dim, dim), dtype=complex)

    for i, bi in enumerate(basis):
        for j, bj in enumerate(basis):
            Sji = gaussian_overlap_general(
                bj.q, bj.p, bj.A,
                bi.q, bi.p, bi.A,
            )
            rho += (
                C[i]
                * np.conj(C[j])
                * Sji
                * np.outer(vectors[i], np.conj(vectors[j]))
            )

    rho = 0.5*(rho+rho.conj().T)
    tr = float(np.real(np.trace(rho)))

    if normalize:
        if tr <= 0.0:
            raise ValueError("electronic density has non-positive trace.")
        rho = rho/tr

    return rho


def reduced_electronic_density_analytic_ci_diabatic(
    coefficients,
    basis,
    normalize=True,
):
    """Reduced density in the global diabatic basis of the analytic 2D CI model."""
    from .ci2d import analytic_adiabatic_vectors

    vectors = [
        analytic_adiabatic_vectors(b.q)[:, int(b.state)]
        for b in basis
    ]
    return reduced_electronic_density_from_vectors(
        coefficients,
        basis,
        vectors,
        normalize=normalize,
    )


def reduced_electronic_density_graph(
    coefficients,
    basis,
    registry,
    reference_node,
    normalize=True,
):
    """Reduced electronic density matrix in one common graph reference frame."""
    dim = registry.graph.dimension
    vectors = []

    for b in basis:
        local = _one_hot(b.state, dim)
        vectors.append(
            registry.transport_coefficients(
                b.node,
                reference_node,
                local,
            )
        )

    return reduced_electronic_density_from_vectors(
        coefficients,
        basis,
        vectors,
        normalize=normalize,
    )

def density_matrix_populations(rho):
    rho = np.asarray(rho, dtype=complex)
    if rho.ndim != 2 or rho.shape[0] != rho.shape[1]:
        raise ValueError("rho must be square.")
    return np.real(np.diag(rho))


def density_matrix_purity(rho):
    rho = np.asarray(rho, dtype=complex)
    tr = np.trace(rho)
    if abs(tr) <= 1e-15:
        raise ValueError("rho has zero trace.")
    normalized = rho / tr
    return float(np.real(np.trace(normalized @ normalized)))


def exact_reduced_electronic_density_diabatic(psi, dx, dy):
    """Integrate a two-state exact-grid wavefunction over nuclear coordinates."""
    psi = np.asarray(psi, dtype=complex)
    if psi.ndim != 3 or psi.shape[-1] != 2:
        raise ValueError("psi must have shape (nx,ny,2).")

    flat = psi.reshape(-1, psi.shape[-1])
    rho = flat.T @ np.conj(flat) * dx * dy
    rho = 0.5 * (rho + rho.conj().T)
    return rho


def rotate_density_to_frame(rho_diabatic, frame):
    """If diabatic vector = frame @ local_vector, then rho_local=frame^dag rho frame."""
    rho = np.asarray(rho_diabatic, dtype=complex)
    U = np.asarray(frame, dtype=complex)
    if rho.shape != U.shape or rho.shape[0] != rho.shape[1]:
        raise ValueError("rho and frame must be equal-size square matrices.")
    return U.conj().T @ rho @ U


def exact_reference_frame_density(psi, dx, dy, frame, normalize=True):
    rho_d = exact_reduced_electronic_density_diabatic(psi, dx, dy)
    rho = rotate_density_to_frame(rho_d, frame)
    tr = float(np.real(np.trace(rho)))
    if normalize:
        rho = rho / tr
    return 0.5*(rho+rho.conj().T)


def density_matrix_linear_entropy(rho):
    """Linear entropy 1-Tr(rho^2) for a normalized reduced density matrix."""
    purity = density_matrix_purity(rho)
    return float(1.0 - purity)


def density_matrix_von_neumann_entropy(rho, logarithm_base=np.e):
    """Von Neumann entropy -Tr(rho log rho) for a normalized Hermitian density."""
    rho = np.asarray(rho, dtype=complex)
    tr = np.trace(rho)
    if abs(tr) <= 1e-15:
        raise ValueError("rho has zero trace.")
    rho = 0.5*(rho/tr + (rho/tr).conj().T)
    eig = np.linalg.eigvalsh(rho).real
    eig = np.clip(eig, 0.0, 1.0)
    nz = eig[eig > 1e-15]
    entropy = -np.sum(nz*np.log(nz))
    if logarithm_base != np.e:
        entropy /= np.log(float(logarithm_base))
    return float(entropy)
