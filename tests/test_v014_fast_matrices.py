import numpy as np

from gaussian_dynamics.dynamic_gauge_graph import AnalyticCI2DFrameProvider
from gaussian_dynamics.dynamic_graph_aims import DynamicGraphTBF
from gaussian_dynamics.spinor_complete_lvc_v12 import (
    build_spinor_complete_lvc_matrices,
)
from gaussian_dynamics.fast_lvc_matrices_v14 import (
    build_spinor_complete_lvc_matrices_symmetric,
    hermitian_pair_evaluation_count,
    ordered_pair_evaluation_count,
    pair_evaluation_reduction,
)


def _tbf(uid,state,q,p,A):
    return DynamicGraphTBF(
        uid=uid,state=state,
        q=np.asarray(q,float),
        p=np.asarray(p,float),
        A=np.asarray(A,float),
        node=("t",uid),
    )


def test_symmetric_half_build_matches_full_ordered_pair_builder():
    provider=AnalyticCI2DFrameProvider(
        nuclear_mass_au=9.0
    )
    basis=[
        _tbf(0,1,[-0.6,0.2],[0.5,0.1],[[1.2,0.05],[0.05,0.9]]),
        _tbf(1,0,[0.2,0.7],[-0.1,0.4],[[0.8,-0.03],[-0.03,1.3]]),
        _tbf(2,1,[0.8,-0.4],[-0.4,0.3],[[1.5,0.0],[0.0,0.7]]),
    ]

    S0,H0,N0=build_spinor_complete_lvc_matrices(
        basis,provider
    )
    S1,H1,N1=build_spinor_complete_lvc_matrices_symmetric(
        basis,provider
    )

    assert np.allclose(S1,S0,atol=2e-13)
    assert np.allclose(H1,H0,atol=2e-13)
    assert np.allclose(N1,N0,atol=2e-13)


def test_pair_evaluation_count_is_asymptotically_halved():
    assert ordered_pair_evaluation_count(10)==100
    assert hermitian_pair_evaluation_count(10)==55
    assert np.isclose(pair_evaluation_reduction(10),0.45)
    assert pair_evaluation_reduction(100)>0.49
