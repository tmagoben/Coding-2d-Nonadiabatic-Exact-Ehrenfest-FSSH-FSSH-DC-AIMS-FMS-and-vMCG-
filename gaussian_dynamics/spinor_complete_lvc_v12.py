import numpy as np

from .gaussian_general import (
    gaussian_overlap_general,
    kinetic_matrix_element_general,
    basis_time_matrix_element_general,
)
from .lvc_exact_gaussian import exact_lvc_potential_matrix_element


def build_nuclear_overlap_matrix(basis):
    n=len(basis)
    S=np.zeros((n,n),dtype=complex)
    for i in range(n):
        for j in range(n):
            S[i,j]=gaussian_overlap_general(
                basis[i].q,basis[i].p,basis[i].A,
                basis[j].q,basis[j].p,basis[j].A,
            )
    return S


def build_spinor_complete_lvc_matrices(basis, provider):
    r"""Exact LVC S/H for a complete two-state electronic basis on every Gaussian.

    Ansatz:

        Psi(R,t) = sum_k g_k(R,t) c_k(t),

    where each c_k is a two-component vector in the fixed global diabatic basis.

    Flattening (k,a) gives

        S_(ka,lb) = <g_k|g_l> delta_ab

        H_(ka,lb) = <g_k|T|g_l> delta_ab
                    + <g_k|V_ab(R)|g_l>.

    No electronic derivative couplings or gauge choices appear.
    """
    n=len(basis)
    ns=2
    dim=n*ns
    Sfull=np.zeros((dim,dim),dtype=complex)
    Hfull=np.zeros((dim,dim),dtype=complex)
    Snuc=np.zeros((n,n),dtype=complex)

    point=provider.evaluate(np.asarray(basis[0].q,float))
    M=np.asarray(point.mass_matrix,float)
    params=provider.params

    for i in range(n):
        si=slice(ns*i,ns*(i+1))
        for j in range(n):
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

            Snuc[i,j]=Sij
            Sfull[si,sj]=Sij*np.eye(ns)
            Hfull[si,sj]=Tij*np.eye(ns)+Vij

    return Sfull,Hfull,Snuc


def build_spinor_complete_time_matrix(
    basis,
    qdots,
    pdots,
):
    """Moving-basis T matrix in a fixed global diabatic electronic basis."""
    n=len(basis)
    ns=2
    dim=n*ns

    qdots=np.asarray(qdots,float)
    pdots=np.asarray(pdots,float)
    if qdots.shape!=(n,len(basis[0].q)) or pdots.shape!=qdots.shape:
        raise ValueError("kinematic arrays have incompatible shape.")

    T=np.zeros((dim,dim),dtype=complex)

    for i in range(n):
        si=slice(ns*i,ns*(i+1))
        for j in range(n):
            sj=slice(ns*j,ns*(j+1))
            t= basis_time_matrix_element_general(
                basis[i].q,basis[i].p,basis[i].A,
                basis[j].q,basis[j].p,basis[j].A,
                qdots[j],pdots[j],
            )
            T[si,sj]=t*np.eye(ns)

    return T


def coefficients_matrix(flat_coefficients, n_gaussian):
    C=np.asarray(flat_coefficients,dtype=complex)
    if C.shape!=(2*int(n_gaussian),):
        raise ValueError("flattened coefficient length must equal 2*n_gaussian.")
    return C.reshape(int(n_gaussian),2)


def flatten_coefficients(C):
    C=np.asarray(C,dtype=complex)
    if C.ndim!=2 or C.shape[1]!=2:
        raise ValueError("C must have shape (n_gaussian,2).")
    return C.reshape(-1)


def spinor_complete_reduced_density(flat_coefficients, basis, normalize=True):
    C=coefficients_matrix(flat_coefficients,len(basis))
    S=build_nuclear_overlap_matrix(basis)

    # rho_ab = sum_ij C_ia C_jb^* <g_j|g_i>
    rho=C.T @ S.T @ np.conj(C)
    rho=0.5*(rho+rho.conj().T)

    if normalize:
        tr=np.trace(rho)
        if abs(tr)<1e-15:
            raise ValueError("zero density trace.")
        rho=rho/tr
    return rho


def spinor_complete_generalized_norm(flat_coefficients, Sfull):
    C=np.asarray(flat_coefficients,dtype=complex)
    return float(np.real(np.vdot(C,np.asarray(Sfull)@C)))
