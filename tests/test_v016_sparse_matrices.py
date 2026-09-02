import numpy as np
from scipy import sparse

from gaussian_dynamics.dynamic_gauge_graph import AnalyticCI2DFrameProvider
from gaussian_dynamics.dynamic_graph_aims import DynamicGraphTBF
from gaussian_dynamics.locality_graph_v16 import (
    LocalityGraphSettings,
    PersistentGaussianLocalityGraph,
)
from gaussian_dynamics.pair_cache_v15 import (
    GaussianPairCache,
    build_cached_spinor_lvc_matrices,
    build_cached_spinor_time_matrix,
)
from gaussian_dynamics.sparse_pair_matrices_v16 import (
    build_sparse_spinor_lvc_matrices,
    build_sparse_spinor_time_matrix,
    sparse_metric_compatible_connection,
    sparse_moving_basis_midpoint_cayley_step,
)
from gaussian_dynamics.moving_graph_gaussian import (
    metric_compatible_basis_connection,
)
from gaussian_dynamics.moving_basis_v12 import (
    moving_basis_midpoint_cayley_step,
)


def _basis():
    return [
        DynamicGraphTBF(
            uid=0,state=1,
            q=np.array([-0.5,0.2]),
            p=np.array([0.5,0.1]),
            A=np.array([[1.2,0.05],[0.05,0.9]]),
            node=("a",0),
        ),
        DynamicGraphTBF(
            uid=1,state=0,
            q=np.array([0.2,0.4]),
            p=np.array([-0.1,0.3]),
            A=np.array([[0.9,-0.03],[-0.03,1.3]]),
            node=("b",1),
        ),
        DynamicGraphTBF(
            uid=2,state=1,
            q=np.array([0.7,-0.5]),
            p=np.array([0.2,-0.4]),
            A=np.array([[1.4,0.0],[0.0,0.8]]),
            node=("c",2),
        ),
    ]


def test_sparse_builder_matches_dense_when_all_edges_active():
    basis=_basis()
    provider=AnalyticCI2DFrameProvider(nuclear_mass_au=9.0)

    graph=PersistentGaussianLocalityGraph(
        LocalityGraphSettings(
            enter_overlap=1e-14,
            exit_overlap=5e-15,
        )
    )
    update=graph.update(basis)
    assert len(update.active_edges)==3

    sparse_mats=build_sparse_spinor_lvc_matrices(
        update,provider
    )

    dense_cache=GaussianPairCache(basis)
    S,H,Snuc=build_cached_spinor_lvc_matrices(
        dense_cache,provider
    )

    assert np.allclose(
        sparse_mats.S.toarray(),S,atol=3e-13
    )
    assert np.allclose(
        sparse_mats.H.toarray(),H,atol=3e-13
    )
    assert np.allclose(
        sparse_mats.Snuc.toarray(),Snuc,atol=3e-13
    )


def test_sparse_time_and_cayley_match_dense_full_graph():
    basis0=_basis()
    basis1=_basis()
    for b in basis1:
        b.q=b.q+np.array([1e-3,-5e-4])

    provider=AnalyticCI2DFrameProvider(nuclear_mass_au=9.0)

    graph0=PersistentGaussianLocalityGraph(
        LocalityGraphSettings(
            enter_overlap=1e-14,
            exit_overlap=5e-15,
        )
    )
    graph1=PersistentGaussianLocalityGraph(
        LocalityGraphSettings(
            enter_overlap=1e-14,
            exit_overlap=5e-15,
        )
    )
    u0=graph0.update(basis0)
    u1=graph1.update(basis1)

    m0=build_sparse_spinor_lvc_matrices(u0,provider)
    m1=build_sparse_spinor_lvc_matrices(u1,provider)

    qdots=np.array([
        [0.05,0.01],
        [0.02,-0.03],
        [-0.01,0.04],
    ])
    pdots=np.array([
        [-0.02,0.01],
        [0.01,0.02],
        [0.03,-0.02],
    ])

    # Midpoint basis is sufficiently close that all graph edges remain active.
    mid=[]
    for a,b in zip(basis0,basis1):
        mid.append(
            DynamicGraphTBF(
                uid=a.uid,state=a.state,
                q=0.5*(a.q+b.q),
                p=0.5*(a.p+b.p),
                A=a.A.copy(),
                node=a.node,
            )
        )
    gm=PersistentGaussianLocalityGraph(
        LocalityGraphSettings(
            enter_overlap=1e-14,
            exit_overlap=5e-15,
        )
    )
    um=gm.update(mid)

    Ts=build_sparse_spinor_time_matrix(
        um,qdots,pdots
    )
    Tcompat_s=sparse_metric_compatible_connection(
        m0.S,m1.S,0.001,Ts
    )

    # Dense reference.
    cd0=GaussianPairCache(basis0)
    cd1=GaussianPairCache(basis1)
    cdm=GaussianPairCache(mid)
    S0,H0,_=build_cached_spinor_lvc_matrices(cd0,provider)
    S1,H1,_=build_cached_spinor_lvc_matrices(cd1,provider)
    Td=build_cached_spinor_time_matrix(cdm,qdots,pdots)
    Tcompat_d=metric_compatible_basis_connection(
        S0,S1,0.001,seed=Td
    )

    assert np.allclose(
        Tcompat_s.toarray(),Tcompat_d,atol=5e-12
    )

    rng=np.random.default_rng(7)
    C=rng.normal(size=6)+1j*rng.normal(size=6)

    Cs=sparse_moving_basis_midpoint_cayley_step(
        C,m0.S,m0.H,m1.S,m1.H,Tcompat_s,0.001
    )
    Cd=moving_basis_midpoint_cayley_step(
        C,S0,H0,S1,H1,Tcompat_d,0.001
    )
    assert np.allclose(Cs,Cd,atol=2e-11)
