import numpy as np

from gaussian_dynamics.dynamic_graph_aims import DynamicGraphTBF
from gaussian_dynamics.dynamic_gauge_graph import AnalyticCI2DFrameProvider
from gaussian_dynamics.spinor_complete_lvc_v12 import (
    build_spinor_complete_lvc_matrices,
    spinor_complete_reduced_density,
    spinor_complete_generalized_norm,
)
from gaussian_dynamics.spinor_complete_dynamics_v12 import (
    initialize_spinor_complete_coefficients,
    run_spinor_complete_lvc_gaussians,
)


def _tbf(uid,state,q,p,A):
    return DynamicGraphTBF(
        uid=uid,
        state=state,
        q=np.asarray(q,float),
        p=np.asarray(p,float),
        A=np.asarray(A,float),
        node=("n",uid),
    )


def test_spinor_complete_blocks_are_hermitian():
    provider=AnalyticCI2DFrameProvider(nuclear_mass_au=8.0)
    basis=[
        _tbf(0,1,[-0.5,0.4],[0.7,0.2],1.2*np.eye(2)),
        _tbf(1,0,[0.4,0.7],[-0.1,0.5],np.diag([0.7,1.4])),
    ]

    S,H,Snuc=build_spinor_complete_lvc_matrices(basis,provider)

    assert S.shape==(4,4)
    assert H.shape==(4,4)
    assert Snuc.shape==(2,2)
    assert np.allclose(S,S.conj().T,atol=1e-12)
    assert np.allclose(H,H.conj().T,atol=1e-12)


def test_center_adiabatic_initialization_has_unit_norm_for_one_gaussian():
    provider=AnalyticCI2DFrameProvider(nuclear_mass_au=20.0)
    basis=[
        _tbf(0,1,[0.6,0.4],[0.2,0.3],1.1*np.eye(2))
    ]
    C=initialize_spinor_complete_coefficients(basis,provider)
    S,_,_=build_spinor_complete_lvc_matrices(basis,provider)

    assert np.isclose(spinor_complete_generalized_norm(C,S),1.0)

    rho=spinor_complete_reduced_density(C,basis,normalize=True)
    assert np.isclose(np.trace(rho),1.0)
    assert np.isclose(np.linalg.matrix_rank(rho,tol=1e-10),1)


def test_spinor_complete_short_run_conserves_norm():
    provider=AnalyticCI2DFrameProvider(nuclear_mass_au=20.0)
    basis=[
        _tbf(0,1,[0.55,0.45],[0.6,0.8],1.2*np.eye(2))
    ]

    out=run_spinor_complete_lvc_gaussians(
        basis,
        provider=provider,
        dt=2e-4,
        steps=20,
        spawn_action_threshold=1e9,
        max_basis=1,
        store_every=1,
    )

    drift=max(abs(r["norm"]-1.0) for r in out["records"])
    assert drift < 2e-6

    rho=spinor_complete_reduced_density(
        out["final_coefficients"],
        out["final_basis"],
        normalize=True,
    )
    assert np.isclose(np.trace(rho),1.0)
