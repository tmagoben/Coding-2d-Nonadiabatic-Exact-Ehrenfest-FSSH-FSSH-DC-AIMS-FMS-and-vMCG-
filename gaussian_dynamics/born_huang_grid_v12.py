from dataclasses import dataclass
import numpy as np

from .ci2d import LVC2DParameters, diabatic_potential_2d
from .gaussian_nd import gaussian_nd, gaussian_nd_gradient, gaussian_nd_time_derivative
from .gaussian_general import (
    gaussian_overlap_general,
    basis_time_matrix_element_general,
)


@dataclass(frozen=True)
class BornHuangGrid2D:
    x: np.ndarray
    X: np.ndarray
    Y: np.ndarray
    points: np.ndarray
    dx: float
    frame: np.ndarray
    derivative_frame: np.ndarray
    diabatic_potential: np.ndarray
    mass: float
    params: LVC2DParameters

    @property
    def area(self):
        return self.dx*self.dx

    @property
    def grid_n(self):
        return len(self.x)


def build_born_huang_grid_2d(
    grid_n=64,
    half_width=4.0,
    mass=5.0,
    params=LVC2DParameters(),
):
    """Precompute the analytic two-state adiabatic frame and its derivatives.

    The half-grid shift is the same convention used by the exact 2D benchmark, so the
    exact CI at (0,0) is not a grid point for even `grid_n`.
    """
    n=int(grid_n)
    if n<=0:
        raise ValueError("grid_n must be positive.")

    L=float(half_width)
    dx=2.0*L/n
    x=-L+(np.arange(n)+0.5)*dx
    X,Y=np.meshgrid(x,x,indexing="ij")
    points=np.stack([X,Y],axis=-1)

    denom=(params.kappa*X)**2+(params.lam*Y)**2
    if np.any(denom<=1e-28):
        raise ValueError(
            "Grid contains the exact conical intersection; use the half-shifted "
            "even grid convention."
        )

    theta=np.arctan2(params.lam*Y,params.kappa*X)
    s=np.sin(0.5*theta)
    c=np.cos(0.5*theta)

    frame=np.zeros((n,n,2,2),dtype=complex)
    frame[...,0,0]=-s
    frame[...,1,0]= c
    frame[...,0,1]= c
    frame[...,1,1]= s

    # d_01 = 1/2 kappa lambda (-y,x) / denom
    a=np.zeros((n,n,2),dtype=float)
    pref=0.5*params.kappa*params.lam/denom
    a[...,0]=-pref*Y
    a[...,1]= pref*X

    # D_alpha[i,j] = <phi_i|partial_alpha phi_j>
    D=np.zeros((n,n,2,2,2),dtype=float)
    D[...,0,1,:]=a
    D[...,1,0,:]=-a

    # derivative_frame[..., alpha, diabatic_component, adiabatic_state]
    derivative_frame=np.zeros((n,n,2,2,2),dtype=complex)
    for alpha in range(2):
        derivative_frame[...,alpha,:,:]=np.einsum(
            "...ac,...cb->...ab",
            frame,
            D[...,alpha],
        )

    V=diabatic_potential_2d(X,Y,params)

    return BornHuangGrid2D(
        x=x,
        X=X,
        Y=Y,
        points=points,
        dx=dx,
        frame=frame,
        derivative_frame=derivative_frame,
        diabatic_potential=V,
        mass=float(mass),
        params=params,
    )


def born_huang_basis_wavefunctions(basis, grid):
    """Return chi_i(R)=g_i(R) Phi_ai(R) without computing nuclear gradients."""
    n=len(basis)
    ng=grid.grid_n
    chi=np.zeros((n,ng,ng,2),dtype=complex)

    for i,b in enumerate(basis):
        state=int(b.state)
        g=gaussian_nd(
            grid.points,
            b.q,b.p,b.A,
        )
        phi=grid.frame[..., :, state]
        chi[i]=g[...,None]*phi

    return chi


def born_huang_basis_fields(basis, grid):
    """Return chi_i(R) and grad chi_i(R) in the diabatic electronic basis."""
    n=len(basis)
    ng=grid.grid_n
    chi=born_huang_basis_wavefunctions(basis,grid)
    grad=np.zeros((n,ng,ng,2,2),dtype=complex)

    for i,b in enumerate(basis):
        state=int(b.state)
        g=gaussian_nd(
            grid.points,
            b.q,b.p,b.A,
        )
        gg=gaussian_nd_gradient(
            grid.points,
            b.q,b.p,b.A,
        )
        phi=grid.frame[..., :, state]
        dphi=grid.derivative_frame[..., :, :, state]

        for alpha in range(2):
            grad[i,...,alpha,:]=(
                gg[...,alpha,None]*phi
                +g[...,None]*dphi[...,alpha,:]
            )

    return chi,grad


def apply_spectral_kinetic_to_basis_fields(chi, grid):
    """Apply the same periodic spectral kinetic operator used by the exact benchmark."""
    chi=np.asarray(chi,dtype=complex)
    if chi.ndim!=4 or chi.shape[-1]!=2:
        raise ValueError("chi must have shape (nbasis,nx,ny,2).")

    nx,ny=chi.shape[1:3]
    kx=2.0*np.pi*np.fft.fftfreq(nx,d=grid.dx)
    ky=2.0*np.pi*np.fft.fftfreq(ny,d=grid.dx)
    KX,KY=np.meshgrid(kx,ky,indexing="ij")
    kinetic_factor=0.5*(KX**2+KY**2)/grid.mass

    chi_k=np.fft.fftn(chi,axes=(1,2))
    return np.fft.ifftn(
        kinetic_factor[None,...,None]*chi_k,
        axes=(1,2),
    )


def build_born_huang_matrices(basis, grid):
    r"""Projected S/H for chi_i(R)=g_i(R) Phi_ai(R) on the exact 2D grid.

    This version intentionally applies the **global diabatic Hamiltonian directly** to
    the coordinate-dependent basis fields.  The kinetic operator is the same periodic
    FFT operator used by `exact2d.py`.

    Consequences:
    - first- and second-derivative electronic couplings do not need to be inserted
      separately;
    - the CI singularity is never evaluated as an explicit second-derivative coupling;
    - branch-cut/geometric-phase structure present in the basis field is handled by
      the same discrete kinetic operator as the exact reference.

    The resulting projected Hamiltonian is exact within:
      1. the selected finite Gaussian basis;
      2. the selected finite periodic 2D grid.
    """
    n=len(basis)
    chi=born_huang_basis_wavefunctions(basis,grid)

    chi_flat=chi.reshape(n,-1)
    S=(chi_flat.conj()@chi_flat.T)*grid.area
    S=0.5*(S+S.conj().T)

    Tchi=apply_spectral_kinetic_to_basis_fields(chi,grid)
    Vchi=np.einsum(
        "...ab,j...b->j...a",
        grid.diabatic_potential,
        chi,
    )
    Hchi=Tchi+Vchi

    H=(chi_flat.conj()@Hchi.reshape(n,-1).T)*grid.area
    H=0.5*(H+H.conj().T)

    return S,H


def born_huang_basis_time_matrix_grid(
    basis,
    grid,
    qdots,
    pdots,
):
    r"""Projected moving-basis matrix on the same grid as S/H.

    Since Phi_a(R) depends on R but not on the moving TBF center parameters,

        dot chi_j(R,t) = dot g_j(R,t) Phi_aj(R).
    """
    n=len(basis)
    qdots=np.asarray(qdots,float)
    pdots=np.asarray(pdots,float)

    if qdots.shape!=(n,2) or pdots.shape!=(n,2):
        raise ValueError("qdots/pdots must have shape (nbasis,2).")

    chi=born_huang_basis_wavefunctions(basis,grid)
    dotchi=np.zeros_like(chi)

    for j,b in enumerate(basis):
        gdot=gaussian_nd_time_derivative(
            grid.points,
            b.q,b.p,b.A,
            qdots[j],pdots[j],
        )
        phi=grid.frame[..., :, int(b.state)]
        dotchi[j]=gdot[...,None]*phi

    return (
        chi.reshape(n,-1).conj()
        @dotchi.reshape(n,-1).T
    )*grid.area

def born_huang_basis_time_matrix(
    basis,
    qdots,
    pdots,
):
    r"""Return <g_i Phi_ai | d/dt(g_j Phi_aj)>.

    Phi_aj(R) depends on the integration coordinate R but not on the moving Gaussian
    center parameters q_j(t), p_j(t).  Therefore

        T_ij = delta_ai,aj <g_i|dot g_j>.
    """
    n=len(basis)
    qdots=np.asarray(qdots,float)
    pdots=np.asarray(pdots,float)

    T=np.zeros((n,n),dtype=complex)

    for i in range(n):
        for j in range(n):
            if int(basis[i].state)!=int(basis[j].state):
                continue
            T[i,j]=basis_time_matrix_element_general(
                basis[i].q,basis[i].p,basis[i].A,
                basis[j].q,basis[j].p,basis[j].A,
                qdots[j],pdots[j],
            )
    return T


def reconstruct_born_huang_wavefunction(coefficients,basis,grid):
    C=np.asarray(coefficients,dtype=complex)
    if C.shape!=(len(basis),):
        raise ValueError("coefficient vector length mismatch.")

    chi=born_huang_basis_wavefunctions(basis,grid)
    return np.einsum("i,i...a->...a",C,chi)


def born_huang_reduced_density(coefficients,basis,grid,normalize=True):
    psi=reconstruct_born_huang_wavefunction(
        coefficients,basis,grid
    )
    flat=psi.reshape(-1,2)
    rho=(flat.T @ np.conj(flat))*grid.area
    rho=0.5*(rho+rho.conj().T)

    if normalize:
        tr=np.trace(rho)
        if abs(tr)<1e-15:
            raise ValueError("zero reduced-density trace.")
        rho=rho/tr
    return rho
