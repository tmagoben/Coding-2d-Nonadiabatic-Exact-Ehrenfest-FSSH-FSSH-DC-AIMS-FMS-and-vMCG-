import numpy as np

from gaussian_dynamics.ci2d import diabatic_potential_2d
from gaussian_dynamics.dynamic_graph_aims import DynamicGraphTBF
from gaussian_dynamics.dynamic_gauge_graph import AnalyticCI2DFrameProvider
from gaussian_dynamics.gaussian_nd import gaussian_nd, gaussian_nd_gradient
from gaussian_dynamics.lvc_exact_gaussian import (
    center_adiabatic_spinor,
    center_spinor_time_derivative,
    exact_lvc_pair_result,
    build_exact_lvc_gaussian_matrices,
    exact_lvc_basis_time_matrix,
)


def _grid(n=190, L=7.0):
    dx = 2*L/n
    x = -L + (np.arange(n)+0.5)*dx
    X,Y = np.meshgrid(x,x,indexing="ij")
    R = np.stack([X,Y],axis=-1)
    return X,Y,R,dx


def _tbf(uid,state,q,p,A):
    return DynamicGraphTBF(
        uid=uid,
        state=state,
        q=np.asarray(q,float),
        p=np.asarray(p,float),
        A=np.asarray(A,float),
        node=("unused",uid),
    )


def test_exact_lvc_pair_matches_direct_2d_quadrature():
    provider=AnalyticCI2DFrameProvider(nuclear_mass_au=7.0)

    bi=_tbf(
        0,1,
        [-0.55,0.35],
        [0.7,-0.2],
        [[1.25,0.08],[0.08,0.85]],
    )
    bj=_tbf(
        1,0,
        [0.45,-0.25],
        [-0.15,0.55],
        [[0.75,-0.05],[-0.05,1.15]],
    )

    result=exact_lvc_pair_result(bi,bj,provider)

    X,Y,R,dx=_grid()
    gi=gaussian_nd(R,bi.q,bi.p,bi.A)
    gj=gaussian_nd(R,bj.q,bj.p,bj.A)

    ui=center_adiabatic_spinor(bi,provider)
    uj=center_adiabatic_spinor(bj,provider)

    V=diabatic_potential_2d(X,Y,provider.params)
    v_scalar=np.einsum(
        "a,...ab,b->...",
        np.conj(ui),
        V,
        uj,
    )
    V_num=np.vdot(gi,v_scalar*gj)*dx*dx

    grad_i=gaussian_nd_gradient(R,bi.q,bi.p,bi.A)
    grad_j=gaussian_nd_gradient(R,bj.q,bj.p,bj.A)
    T_nuc=0j
    Minv=np.eye(2)/provider.nuclear_mass_au
    for a in range(2):
        for b in range(2):
            T_nuc += (
                0.5*Minv[a,b]
                *np.vdot(grad_i[...,a],grad_j[...,b])
                *dx*dx
            )
    T_num=T_nuc*np.vdot(ui,uj)

    S_num=np.vdot(gi,gj)*np.vdot(ui,uj)*dx*dx

    assert abs(result.overlap-S_num) < 5e-10
    assert abs(result.potential-V_num) < 5e-10
    assert abs(result.kinetic-T_num) < 5e-10
    assert abs(result.total-(T_num+V_num)) < 8e-10


def test_exact_lvc_matrices_are_hermitian_for_mixed_states_and_widths():
    provider=AnalyticCI2DFrameProvider(nuclear_mass_au=10.0)
    basis=[
        _tbf(0,1,[-0.7,0.4],[0.8,0.1],[[1.2,0.1],[0.1,0.9]]),
        _tbf(1,0,[0.2,0.6],[0.1,0.7],[[0.7,-0.05],[-0.05,1.4]]),
        _tbf(2,1,[0.8,-0.3],[-0.5,0.2],[[1.5,0.0],[0.0,0.6]]),
    ]

    S,H=build_exact_lvc_gaussian_matrices(basis,provider)

    assert np.allclose(S,S.conj().T,atol=2e-12)
    assert np.allclose(H,H.conj().T,atol=2e-12)
    assert np.min(np.linalg.eigvalsh(S)) > 0.0


def test_center_spinor_time_derivative_matches_finite_difference():
    provider=AnalyticCI2DFrameProvider(nuclear_mass_au=12.0)
    b=_tbf(0,1,[0.8,0.65],[0.4,-0.2],1.1*np.eye(2))
    qdot=np.array([0.17,-0.11])

    analytic=center_spinor_time_derivative(b,provider,qdot=qdot)

    h=1e-6
    bp=_tbf(0,1,b.q+h*qdot,b.p,b.A)
    bm=_tbf(0,1,b.q-h*qdot,b.p,b.A)
    fd=(
        center_adiabatic_spinor(bp,provider)
        -center_adiabatic_spinor(bm,provider)
    )/(2*h)

    assert np.allclose(analytic,fd,atol=2e-9)


def test_exact_basis_connection_satisfies_continuous_metric_identity():
    provider=AnalyticCI2DFrameProvider(nuclear_mass_au=9.0)
    basis=[
        _tbf(0,1,[0.8,0.55],[0.4,0.2],[[1.2,0.06],[0.06,0.9]]),
        _tbf(1,0,[0.25,0.85],[-0.1,0.35],[[0.8,-0.04],[-0.04,1.3]]),
    ]

    qdots=np.array([[0.07,-0.03],[-0.02,0.05]])
    pdots=np.array([[0.01,0.02],[-0.03,0.01]])

    T=exact_lvc_basis_time_matrix(
        basis,provider,qdots,pdots
    )

    h=2e-6

    def shifted(sign):
        out=[]
        for k,b in enumerate(basis):
            out.append(_tbf(
                b.uid,b.state,
                b.q+sign*h*qdots[k],
                b.p+sign*h*pdots[k],
                b.A,
            ))
        return out

    Sp,_=build_exact_lvc_gaussian_matrices(shifted(+1),provider)
    Sm,_=build_exact_lvc_gaussian_matrices(shifted(-1),provider)
    Sdot=(Sp-Sm)/(2*h)

    assert np.allclose(Sdot,T+T.conj().T,atol=2e-8)
