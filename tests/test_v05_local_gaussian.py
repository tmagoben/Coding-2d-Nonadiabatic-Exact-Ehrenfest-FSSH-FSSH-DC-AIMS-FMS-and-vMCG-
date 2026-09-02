import numpy as np

from gaussian_dynamics.gaussian_nd import gaussian_nd, gaussian_nd_gradient
from gaussian_dynamics.local_gaussian_nd import (
    overlap_centroid_equal_width,
    gradient_matrix_element_equal_width,
    kinetic_matrix_element_equal_width,
    LocalAdiabaticTBF,
    local_matrices,
    basis_time_matrix_element_equal_width,
)
from gaussian_dynamics.benchmark_provider_nd import LVC2DGeneralizedProvider


def grid(n=180,L=7.0):
    dx=2*L/n
    x=-L+(np.arange(n)+0.5)*dx
    X,Y=np.meshgrid(x,x,indexing="ij")
    P=np.stack([X,Y],axis=-1)
    return P,dx


def test_gradient_and_kinetic_matrix_elements_match_2d_quadrature():
    P,dx=grid()
    A=np.array([[1.1,0.12],[0.12,0.9]])
    qi=np.array([-0.7,0.3]); pi=np.array([0.6,-0.2])
    qj=np.array([0.5,-0.4]); pj=np.array([-0.1,0.7])
    M=np.array([[3.0,0.25],[0.25,2.0]])

    gi=gaussian_nd(P,qi,pi,A)
    gj=gaussian_nd(P,qj,pj,A)
    grad_i=gaussian_nd_gradient(P,qi,pi,A)
    grad_j=gaussian_nd_gradient(P,qj,pj,A)

    numeric_G=np.array([
        np.vdot(gi,grad_j[...,k])*dx*dx
        for k in range(2)
    ])
    analytic_G=gradient_matrix_element_equal_width(qi,pi,qj,pj,A)

    B=np.linalg.inv(M)
    numeric_T=0.0+0.0j
    for a in range(2):
        for b in range(2):
            numeric_T += 0.5*B[a,b]*np.vdot(grad_i[...,a],grad_j[...,b])*dx*dx

    analytic_T=kinetic_matrix_element_equal_width(qi,pi,qj,pj,A,M)

    assert np.allclose(numeric_G,analytic_G,atol=2e-10)
    assert abs(numeric_T-analytic_T) < 2e-10


def test_local_matrices_are_hermitian():
    provider=LVC2DGeneralizedProvider(nuclear_mass_au=20.0)
    A=np.eye(2)

    basis=[
        LocalAdiabaticTBF(0,np.array([-1.1,0.7]),np.array([0.4,0.1]),A),
        LocalAdiabaticTBF(1,np.array([0.8,0.9]),np.array([-0.2,0.3]),A),
    ]

    S,H,T=local_matrices(basis,provider)

    assert np.allclose(S,S.conj().T,atol=1e-12)
    assert np.allclose(H,H.conj().T,atol=1e-11)


def test_basis_time_overlap_identity_by_finite_difference():
    provider=LVC2DGeneralizedProvider(nuclear_mass_au=20.0)
    A=np.eye(2)

    basis=[
        LocalAdiabaticTBF(1,np.array([-1.0,0.8]),np.array([0.5,0.1]),A),
        LocalAdiabaticTBF(1,np.array([-0.2,0.9]),np.array([0.2,-0.1]),A),
    ]

    S,H,T=local_matrices(basis,provider)

    from gaussian_dynamics.local_gaussian_nd import tbf_guidance, local_overlap_element

    h=1e-5
    plus=[]; minus=[]
    for b in basis:
        qd,pd=tbf_guidance(b,provider)
        plus.append(LocalAdiabaticTBF(b.state,b.q+h*qd,b.p+h*pd,b.A))
        minus.append(LocalAdiabaticTBF(b.state,b.q-h*qd,b.p-h*pd,b.A))

    Sp=np.array([[local_overlap_element(a,b) for b in plus] for a in plus])
    Sm=np.array([[local_overlap_element(a,b) for b in minus] for a in minus])
    Sdot=(Sp-Sm)/(2*h)

    assert np.allclose(Sdot,T+T.conj().T,atol=2e-9)
