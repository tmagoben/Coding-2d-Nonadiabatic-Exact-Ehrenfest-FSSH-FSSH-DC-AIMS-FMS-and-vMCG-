from dataclasses import dataclass
import numpy as np

from .dynamic_graph_aims import DynamicGraphTBF
from .gaussian_nd import gaussian_nd
from .spinor_complete_lvc_v12 import (
    build_spinor_complete_lvc_matrices,
    flatten_coefficients,
)


@dataclass(frozen=True)
class InitialProjectionResult:
    coefficients: np.ndarray
    projected_wavefunction: np.ndarray
    target_norm: float
    projected_norm: float
    residual_norm: float
    relative_residual: float
    fidelity: float
    condition_number: float


def make_shifted_initial_gaussian_bank(
    q0,
    p0,
    A,
    state,
    shifts=((0.0,0.0),(0.18,0.0),(-0.18,0.0),(0.0,0.18),(0.0,-0.18)),
    uid_start=0,
):
    q0=np.asarray(q0,float)
    p0=np.asarray(p0,float)
    A=np.asarray(A,float)

    basis=[]
    for k,shift in enumerate(shifts):
        shift=np.asarray(shift,float)
        if shift.shape!=q0.shape:
            raise ValueError("every shift must match q0 dimension.")
        basis.append(
            DynamicGraphTBF(
                uid=int(uid_start+k),
                state=int(state),
                q=q0+shift,
                p=p0.copy(),
                A=A.copy(),
                node=("initial_bank",int(uid_start+k)),
            )
        )
    return basis


def project_grid_wavefunction_to_spinor_complete_basis(
    psi_target,
    points,
    dx,
    basis,
    provider,
):
    r"""Least-squares Hilbert-space projection onto g_i(R)|a_d>.

    The normal equations are

        S C = b,

    with

        b_(i,a) = <g_i a_d | Psi_target>.

    For the spinor-complete basis S = S_nuclear \otimes I_2.
    """
    psi=np.asarray(psi_target,dtype=complex)
    points=np.asarray(points,dtype=float)
    dx=float(dx)

    if psi.shape!=points.shape[:-1]+(2,):
        raise ValueError("psi_target and points have incompatible shapes.")

    Sfull,_,_=build_spinor_complete_lvc_matrices(
        basis,provider
    )

    n=len(basis)
    b=np.zeros((n,2),dtype=complex)

    for i,tbf in enumerate(basis):
        g=gaussian_nd(points,tbf.q,tbf.p,tbf.A)
        for a in range(2):
            b[i,a]=np.vdot(g,psi[...,a])*dx*dx

    rhs=flatten_coefficients(b)

    # Solve the projection equations.  If the bank is deliberately redundant, use
    # least squares rather than an explicit inverse.
    C,_,_,_=np.linalg.lstsq(Sfull,rhs,rcond=1e-12)

    projected=np.zeros_like(psi)
    Cmat=C.reshape(n,2)
    for i,tbf in enumerate(basis):
        g=gaussian_nd(points,tbf.q,tbf.p,tbf.A)
        projected += g[...,None]*Cmat[i][None,None,:]

    target_norm=float(np.sum(np.abs(psi)**2)*dx*dx)
    projected_norm=float(np.sum(np.abs(projected)**2)*dx*dx)
    residual=float(np.sum(np.abs(projected-psi)**2)*dx*dx)

    overlap=np.sum(np.conj(psi)*projected)*dx*dx
    fidelity=float(
        abs(overlap)**2
        /max(target_norm*projected_norm,1e-30)
    )

    return InitialProjectionResult(
        coefficients=C,
        projected_wavefunction=projected,
        target_norm=target_norm,
        projected_norm=projected_norm,
        residual_norm=residual,
        relative_residual=residual/max(target_norm,1e-30),
        fidelity=fidelity,
        condition_number=float(np.linalg.cond(Sfull)),
    )
