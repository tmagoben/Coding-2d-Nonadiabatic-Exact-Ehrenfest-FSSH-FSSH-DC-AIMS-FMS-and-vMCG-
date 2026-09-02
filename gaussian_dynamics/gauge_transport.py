import numpy as np


def align_subspace(reference, current):
    """Unitary/orthogonal Procrustes alignment of current to reference.

    Both arrays have shape (ambient_dimension, subspace_dimension) and orthonormal
    columns.  The returned current @ Q minimizes the Frobenius distance to reference.
    """
    reference = np.asarray(reference, dtype=complex)
    current = np.asarray(current, dtype=complex)

    if reference.shape != current.shape:
        raise ValueError("reference and current must have equal shape.")

    M = current.conj().T @ reference
    U, _, Vh = np.linalg.svd(M)
    Q = U @ Vh
    return current @ Q, Q


def subspace_projector(U):
    U = np.asarray(U, dtype=complex)
    return U @ U.conj().T


def projector_distance(U, V):
    P = subspace_projector(U)
    Q = subspace_projector(V)
    return float(np.linalg.norm(P-Q, ord="fro"))
