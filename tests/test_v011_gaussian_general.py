import numpy as np

from gaussian_dynamics.gaussian_nd import gaussian_nd, gaussian_nd_gradient
from gaussian_dynamics.gaussian_general import (
    gaussian_overlap_general,
    gaussian_cross_centroid,
    real_overlap_saddle_point,
    gradient_matrix_element_general,
    kinetic_matrix_element_general,
    basis_time_matrix_element_general,
)


def grid(n=190,L=7.0):
    dx=2*L/n
    x=-L+(np.arange(n)+0.5)*dx
    X,Y=np.meshgrid(x,x,indexing="ij")
    P=np.stack([X,Y],axis=-1)
    return P,dx


def test_unequal_width_overlap_gradient_and_kinetic_match_quadrature():
    P,dx=grid()
    Ai=np.array([[1.25,0.10],[0.10,0.80]])
    Aj=np.array([[0.75,-0.08],[-0.08,1.10]])
    qi=np.array([-0.6,0.3]); pi=np.array([0.7,-0.15])
    qj=np.array([0.5,-0.45]); pj=np.array([-0.2,0.65])
    M=np.array([[2.7,0.2],[0.2,1.9]])

    gi=gaussian_nd(P,qi,pi,Ai)
    gj=gaussian_nd(P,qj,pj,Aj)
    grad_i=gaussian_nd_gradient(P,qi,pi,Ai)
    grad_j=gaussian_nd_gradient(P,qj,pj,Aj)

    S_num=np.vdot(gi,gj)*dx*dx
    S=gaussian_overlap_general(qi,pi,Ai,qj,pj,Aj)
    assert abs(S-S_num) < 3e-10

    G_num=np.array([
        np.vdot(gi,grad_j[...,k])*dx*dx
        for k in range(2)
    ])
    G=gradient_matrix_element_general(qi,pi,Ai,qj,pj,Aj)
    assert np.allclose(G,G_num,atol=3e-10)

    Minv=np.linalg.inv(M)
    T_num=0j
    for a in range(2):
        for b in range(2):
            T_num += (
                0.5*Minv[a,b]
                *np.vdot(grad_i[...,a],grad_j[...,b])*dx*dx
            )

    T=kinetic_matrix_element_general(qi,pi,Ai,qj,pj,Aj,M)
    assert abs(T-T_num) < 4e-10


def test_general_formula_reduces_to_equal_width_centroid_and_saddle():
    A=np.array([[1.1,0.1],[0.1,0.9]])
    qi=np.array([-0.4,0.2]); qj=np.array([0.8,-0.1])
    pi=np.array([0.3,0.7]); pj=np.array([-0.2,0.1])

    qc=real_overlap_saddle_point(qi,A,qj,A)
    mu=gaussian_cross_centroid(qi,pi,A,qj,pj,A)

    assert np.allclose(qc,0.5*(qi+qj),atol=1e-13)
    expected=0.5*(qi+qj)+0.5j*np.linalg.solve(A,pj-pi)
    assert np.allclose(mu,expected,atol=1e-13)


def test_basis_time_element_with_width_motion_matches_finite_difference():
    P,dx=grid(n=150,L=6.5)

    Ai=np.array([[1.1,0.08],[0.08,0.85]])
    Aj=np.array([[0.8,-0.04],[-0.04,1.2]])

    qi=np.array([-0.2,0.5]); pi=np.array([0.4,0.1])
    qj=np.array([0.55,-0.3]); pj=np.array([-0.15,0.45])

    qdot=np.array([0.2,-0.1])
    pdot=np.array([0.03,0.04])
    Adot=np.array([[0.05,0.01],[0.01,-0.03]])

    gi=gaussian_nd(P,qi,pi,Ai)

    h=1e-6
    gjp=gaussian_nd(
        P,
        qj+h*qdot,
        pj+h*pdot,
        Aj+h*Adot,
    )
    gjm=gaussian_nd(
        P,
        qj-h*qdot,
        pj-h*pdot,
        Aj-h*Adot,
    )
    gdot_fd=(gjp-gjm)/(2*h)

    num=np.vdot(gi,gdot_fd)*dx*dx
    analytic=basis_time_matrix_element_general(
        qi,pi,Ai,
        qj,pj,Aj,
        qdot,pdot,Adot,
    )

    assert abs(num-analytic) < 2e-8
