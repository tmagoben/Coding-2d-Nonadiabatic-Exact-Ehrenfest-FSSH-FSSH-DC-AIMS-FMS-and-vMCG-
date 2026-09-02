import numpy as np

from gaussian_dynamics.gaussian_nd import (
    gaussian_nd,
    analytic_overlap_equal_width,
    gaussian_nd_laplacian,
)


def make_grid(n=160,L=7.0):
    dx=2*L/n
    x=-L+(np.arange(n)+0.5)*dx
    y=x.copy()
    X,Y=np.meshgrid(x,y,indexing="ij")
    P=np.stack([X,Y],axis=-1)
    return x,y,X,Y,P,dx


def test_multidimensional_gaussian_norm_and_covariance():
    x,y,X,Y,P,dx=make_grid()
    A=np.array([[1.2,0.2],[0.2,0.9]])
    q=np.array([-0.4,0.7])
    p=np.array([0.8,-0.3])
    g=gaussian_nd(P,q,p,A)

    w=dx*dx
    norm=np.sum(np.abs(g)**2)*w
    assert abs(norm-1.0) < 2e-11

    xi=np.stack([X-q[0],Y-q[1]],axis=-1)
    flat=xi.reshape(-1,2)
    dens=np.abs(g).reshape(-1)**2
    cov=flat.T @ (flat*dens[:,None]) * w
    assert np.allclose(cov,0.5*np.linalg.inv(A),atol=2e-10)


def test_multidimensional_analytic_overlap_matches_quadrature():
    x,y,X,Y,P,dx=make_grid()
    A=np.array([[1.1,0.15],[0.15,0.8]])
    qi=np.array([-0.8,0.3]); pi=np.array([0.6,0.2])
    qj=np.array([0.5,-0.4]); pj=np.array([-0.2,0.7])

    gi=gaussian_nd(P,qi,pi,A)
    gj=gaussian_nd(P,qj,pj,A)
    numeric=np.vdot(gi,gj)*dx*dx
    analytic=analytic_overlap_equal_width(qi,pi,qj,pj,A)

    assert abs(numeric-analytic) < 2e-10


def test_analytic_laplacian_matches_fft_for_localized_packet():
    x,y,X,Y,P,dx=make_grid(n=128,L=8.0)
    A=np.array([[1.0,0.12],[0.12,0.85]])
    q=np.array([-0.3,0.4]); p=np.array([0.7,-0.2])
    g=gaussian_nd(P,q,p,A)
    lap=gaussian_nd_laplacian(P,q,p,A)

    k=2*np.pi*np.fft.fftfreq(len(x),d=dx)
    KX,KY=np.meshgrid(k,k,indexing="ij")
    lap_fft=np.fft.ifftn(-(KX**2+KY**2)*np.fft.fftn(g))

    err=np.sqrt(np.sum(np.abs(lap-lap_fft)**2)*dx*dx)
    ref=np.sqrt(np.sum(np.abs(lap)**2)*dx*dx)
    assert err/ref < 2e-9
