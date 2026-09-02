import numpy as np


def _normalized_density(rho):
    rho = np.asarray(rho, dtype=complex)
    if rho.ndim != 2 or rho.shape[0] != rho.shape[1]:
        raise ValueError("rho must be square.")
    rho = 0.5*(rho+rho.conj().T)
    tr = np.trace(rho)
    if abs(tr) < 1e-15:
        raise ValueError("rho has zero trace.")
    return rho/tr


def offdiagonal_coherence(rho):
    rho = _normalized_density(rho)
    if rho.shape != (2,2):
        raise ValueError("offdiagonal_coherence currently expects a two-state density.")
    return complex(rho[0,1])


def coherence_magnitude(rho):
    return float(abs(offdiagonal_coherence(rho)))


def wrapped_phase_difference(a, b):
    """Return arg(a)-arg(b) wrapped to [-pi,pi]."""
    a = complex(a)
    b = complex(b)
    return float(np.angle(a*np.conj(b)))


def coherence_phase_error(candidate_rho, reference_rho, magnitude_floor=1e-8):
    a = offdiagonal_coherence(candidate_rho)
    b = offdiagonal_coherence(reference_rho)

    if abs(a) < magnitude_floor or abs(b) < magnitude_floor:
        return None
    return abs(wrapped_phase_difference(a,b))


def coherence_magnitude_error(candidate_rho, reference_rho):
    return abs(coherence_magnitude(candidate_rho)-coherence_magnitude(reference_rho))


def density_trace_distance(candidate_rho, reference_rho):
    """Trace distance 1/2 ||rho-sigma||_1 for Hermitian density matrices."""
    a = _normalized_density(candidate_rho)
    b = _normalized_density(reference_rho)
    delta = 0.5*((a-b)+(a-b).conj().T)
    eig = np.linalg.eigvalsh(delta)
    return float(0.5*np.sum(np.abs(eig)))


def bloch_vector(rho):
    """Bloch vector (x,y,z) for a normalized 2x2 electronic density matrix."""
    rho = _normalized_density(rho)
    if rho.shape != (2,2):
        raise ValueError("bloch_vector requires a 2x2 density matrix.")

    return np.array([
        2.0*np.real(rho[0,1]),
        -2.0*np.imag(rho[0,1]),
        np.real(rho[0,0]-rho[1,1]),
    ])


def bloch_vector_error(candidate_rho, reference_rho):
    return float(np.linalg.norm(
        bloch_vector(candidate_rho)-bloch_vector(reference_rho)
    ))
