import numpy as np

from .gaussian_general import (
    gaussian_overlap_general,
    kinetic_matrix_element_general,
)
from .lvc_exact_gaussian import (
    exact_lvc_potential_matrix_element,
)


def build_spinor_complete_lvc_matrices_symmetric(
    basis,
    provider,
):
    r"""Hermitian half-build of the exact spinor-complete LVC matrices.

    The v0.12 reference builder evaluates every ordered Gaussian pair (i,j), for N^2
    pair evaluations.

    Hermiticity gives

        S_ji = S_ij^*
        H_ji = H_ij^\dagger.

    v0.14 therefore evaluates only

        i <= j,

    requiring

        N(N+1)/2

    pair evaluations.

    The asymptotic complexity remains O(N^2 d^3), but the leading pair-evaluation
    count approaches one half of the full ordered-pair implementation.
    """
    n=len(basis)
    if n==0:
        raise ValueError("basis cannot be empty.")

    ns=2
    dim=ns*n
    Sfull=np.zeros((dim,dim),dtype=complex)
    Hfull=np.zeros((dim,dim),dtype=complex)
    Snuc=np.zeros((n,n),dtype=complex)

    point=provider.evaluate(
        np.asarray(basis[0].q,float)
    )
    M=np.asarray(point.mass_matrix,float)
    params=provider.params

    eye=np.eye(ns,dtype=complex)

    for i in range(n):
        si=slice(ns*i,ns*(i+1))
        for j in range(i,n):
            sj=slice(ns*j,ns*(j+1))

            Sij=gaussian_overlap_general(
                basis[i].q,basis[i].p,basis[i].A,
                basis[j].q,basis[j].p,basis[j].A,
            )
            Tij=kinetic_matrix_element_general(
                basis[i].q,basis[i].p,basis[i].A,
                basis[j].q,basis[j].p,basis[j].A,
                M,
            )
            Vij=exact_lvc_potential_matrix_element(
                basis[i].q,basis[i].p,basis[i].A,
                basis[j].q,basis[j].p,basis[j].A,
                params=params,
            )

            blockS=Sij*eye
            blockH=Tij*eye+Vij

            Snuc[i,j]=Sij
            Sfull[si,sj]=blockS
            Hfull[si,sj]=blockH

            if i!=j:
                Snuc[j,i]=np.conj(Sij)
                Sfull[sj,si]=blockS.conj().T
                Hfull[sj,si]=blockH.conj().T

    return Sfull,Hfull,Snuc


def hermitian_pair_evaluation_count(n_basis):
    n=int(n_basis)
    if n<0:
        raise ValueError("n_basis cannot be negative.")
    return n*(n+1)//2


def ordered_pair_evaluation_count(n_basis):
    n=int(n_basis)
    if n<0:
        raise ValueError("n_basis cannot be negative.")
    return n*n


def pair_evaluation_reduction(n_basis):
    full=ordered_pair_evaluation_count(n_basis)
    if full==0:
        return 0.0
    half=hermitian_pair_evaluation_count(n_basis)
    return 1.0-half/full
