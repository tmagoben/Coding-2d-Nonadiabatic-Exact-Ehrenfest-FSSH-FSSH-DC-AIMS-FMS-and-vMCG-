import numpy as np

from .gaussian_general import basis_time_matrix_element_general


def nuclear_seed_basis_time_matrix_general(
    basis,
    registry,
    reference_selector,
    qdots,
    pdots,
    Adots=None,
):
    """Nuclear moving-basis seed allowing unequal and time-dependent real widths."""
    n=len(basis)
    qdots=np.asarray(qdots,dtype=float)
    pdots=np.asarray(pdots,dtype=float)

    if qdots.shape != (n,len(basis[0].q)) or pdots.shape != qdots.shape:
        raise ValueError("kinematic arrays have incompatible shape.")

    if Adots is None:
        Adots=[np.zeros_like(b.A) for b in basis]
    if len(Adots)!=n:
        raise ValueError("Adots must contain one matrix per basis function.")

    T=np.zeros((n,n),dtype=complex)

    for i in range(n):
        for j in range(n):
            ref=reference_selector(i,j)
            factors=registry.pair_factors(
                basis[i].node,
                basis[i].electronic_coefficients,
                basis[j].node,
                basis[j].electronic_coefficients,
                ref,
            )
            Tn=basis_time_matrix_element_general(
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
            T[i,j]=Tn*factors["overlap"]

    return T
