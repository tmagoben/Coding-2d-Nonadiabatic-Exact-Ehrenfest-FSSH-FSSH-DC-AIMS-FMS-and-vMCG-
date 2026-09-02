import numpy as np

from .finite_manifold_transport_v233 import (
    CONSUMER_OVERLAP_POLICY_V233,
    certified_transport_from_overlap_v233,
)


def nearest_unitary(overlap):
    """Unitary polar factor of a square overlap matrix.

    If O = U Sigma V^dagger, the nearest unitary to O in Frobenius norm is

        W = U V^dagger.
    """
    result = certified_transport_from_overlap_v233(
        overlap, policy=CONSUMER_OVERLAP_POLICY_V233
    )
    return result.right_to_left_transport.copy()


def current_to_previous_procrustes(overlap):
    """Return Q that rotates the current electronic basis toward the previous basis.

    overlap is

        O_ij = <previous_i | current_j>.

    Transform current states as

        |current'> = |current> Q.

    The Procrustes solution makes O Q positive Hermitian:

        Q = V U^dagger,  for O = U Sigma V^dagger.
    """
    result = certified_transport_from_overlap_v233(
        overlap, policy=CONSUMER_OVERLAP_POLICY_V233
    )
    O = result.overlap
    s = result.singular_values
    Q = result.right_to_left_transport.conj().T
    aligned_overlap = O @ Q

    return Q, aligned_overlap, s


def directional_nac_from_overlap(overlap, displacement):
    """First-order anti-Hermitian derivative coupling along a path segment.

    For a scalar path coordinate s,

        O_ij(s,s+ds)
          = delta_ij + ds * d_ij^(s) + O(ds^2),

    so

        d^(s) ~= (O - O^dagger)/(2 ds).

    This is a first-order local diagnostic, not a replacement for an analytic
    multidimensional NAC vector.
    """
    O = np.asarray(overlap, dtype=complex)
    ds = float(displacement)

    if O.ndim != 2 or O.shape[0] != O.shape[1]:
        raise ValueError("overlap must be square.")
    if ds == 0.0:
        raise ValueError("displacement must be nonzero.")

    return (O - O.conj().T) / (2.0 * ds)


def overlap_unitarity_defect(overlap):
    """||O^dagger O - I||_F, useful for judging finite-step/subspace leakage."""
    O = np.asarray(overlap, dtype=complex)
    I = np.eye(O.shape[1], dtype=complex)
    return float(np.linalg.norm(O.conj().T @ O - I, ord="fro"))


def principal_angles(overlap):
    """Principal angles inferred from singular values of a subspace overlap matrix."""
    s = np.linalg.svd(np.asarray(overlap, dtype=complex), compute_uv=False)
    s = np.clip(s, 0.0, 1.0)
    return np.arccos(s)
