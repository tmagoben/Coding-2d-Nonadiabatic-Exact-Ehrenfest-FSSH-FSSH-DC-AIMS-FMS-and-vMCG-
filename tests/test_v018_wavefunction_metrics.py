import numpy as np

from gaussian_dynamics.wavefunction_metrics_v18 import (
    phase_aligned_fidelity,
    phase_aligned_l2_error,
    nuclear_density_l2_error,
    nuclear_density_total_variation,
    spatial_moments,
    compare_wavefunctions,
)


def _grid(n=24):
    x=np.linspace(-2.0,2.0,n,endpoint=False)
    dx=x[1]-x[0]
    X,Y=np.meshgrid(x,x,indexing="ij")
    points=np.stack([X,Y],axis=-1)
    return points,dx*dx


def _psi(points):
    X=points[...,0]
    Y=points[...,1]
    g=np.exp(-0.8*((X-0.2)**2+(Y+0.1)**2))
    psi=np.zeros(X.shape+(2,),complex)
    psi[...,0]=g
    psi[...,1]=0.3j*g*np.exp(0.2j*X)
    return psi


def test_global_phase_is_removed_exactly():
    points,area=_grid()
    psi=_psi(points)
    shifted=np.exp(1.234j)*psi

    assert abs(
        phase_aligned_fidelity(psi,shifted,area)-1.0
    )<1e-12
    assert phase_aligned_l2_error(
        psi,shifted,area
    )<1e-12
    assert nuclear_density_l2_error(
        psi,shifted,area
    )<1e-12


def test_density_and_moment_metrics_detect_translation():
    points,area=_grid()
    psi=_psi(points)

    shifted_points=points.copy()
    X=points[...,0]
    Y=points[...,1]
    g=np.exp(-0.8*((X-0.5)**2+(Y+0.1)**2))
    shifted=np.zeros_like(psi)
    shifted[...,0]=g
    shifted[...,1]=0.3j*g*np.exp(0.2j*X)

    assert nuclear_density_l2_error(
        psi,shifted,area
    )>0.0
    assert nuclear_density_total_variation(
        psi,shifted,area
    )>0.0

    a=spatial_moments(psi,points,area)
    b=spatial_moments(shifted,points,area)
    assert b["mean"][0]>a["mean"][0]

    out=compare_wavefunctions(
        psi,shifted,points,area
    )
    assert out["fidelity"]<1.0
    assert out["mean_error_l2"]>0.0
