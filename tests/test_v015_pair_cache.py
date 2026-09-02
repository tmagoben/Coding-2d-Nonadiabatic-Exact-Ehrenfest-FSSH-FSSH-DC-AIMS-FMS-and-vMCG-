import numpy as np

from gaussian_dynamics.dynamic_gauge_graph import AnalyticCI2DFrameProvider
from gaussian_dynamics.dynamic_graph_aims import DynamicGraphTBF
from gaussian_dynamics.spinor_complete_lvc_v12 import (
    build_spinor_complete_lvc_matrices,
    build_spinor_complete_time_matrix,
)
from gaussian_dynamics.pair_cache_v15 import (
    GaussianPairCache,
    build_cached_spinor_lvc_matrices,
    build_cached_spinor_time_matrix,
    expand_cached_spinor_lvc_matrices,
    subset_cached_spinor_lvc_matrices,
    v14_factorization_equivalent_for_sh,
    v14_factorization_equivalent_for_time,
    v15_factorization_count_for_snapshot,
)


def _basis():
    return [
        DynamicGraphTBF(
            uid=0,state=1,
            q=np.array([-0.55,0.35]),
            p=np.array([0.7,-0.2]),
            A=np.array([[1.25,0.08],[0.08,0.85]]),
            node=("a",0),
        ),
        DynamicGraphTBF(
            uid=1,state=0,
            q=np.array([0.45,-0.25]),
            p=np.array([-0.15,0.55]),
            A=np.array([[0.75,-0.05],[-0.05,1.15]]),
            node=("b",1),
        ),
        DynamicGraphTBF(
            uid=2,state=1,
            q=np.array([0.15,0.8]),
            p=np.array([0.2,-0.4]),
            A=np.array([[1.4,0.02],[0.02,0.7]]),
            node=("c",2),
        ),
    ]


def test_cached_pair_algebra_matches_reference_sh_and_t():
    basis=_basis()
    provider=AnalyticCI2DFrameProvider(
        nuclear_mass_au=8.0
    )
    qdots=np.array([
        [0.05,0.02],
        [-0.01,0.03],
        [0.02,-0.04],
    ])
    pdots=np.array([
        [-0.02,0.01],
        [0.03,-0.01],
        [-0.01,0.02],
    ])

    S0,H0,N0=build_spinor_complete_lvc_matrices(
        basis,provider
    )
    T0=build_spinor_complete_time_matrix(
        basis,qdots,pdots
    )

    cache=GaussianPairCache(basis)
    S1,H1,N1=build_cached_spinor_lvc_matrices(
        cache,provider
    )
    solves_after_sh=cache.stats.canonical_solves
    T1=build_cached_spinor_time_matrix(
        cache,qdots,pdots
    )

    assert np.allclose(S1,S0,atol=3e-13)
    assert np.allclose(H1,H0,atol=3e-13)
    assert np.allclose(N1,N0,atol=3e-13)
    assert np.allclose(T1,T0,atol=3e-13)

    # T reuses the exact same pair moments: no new dense pair solve is needed.
    assert solves_after_sh==6
    assert cache.stats.canonical_solves==6


def test_reverse_orientation_uses_conjugate_pair_view_without_new_solve():
    basis=_basis()[:2]
    cache=GaussianPairCache(basis)

    a=cache.pair(0,1)
    nsolve=cache.stats.canonical_solves
    b=cache.pair(1,0)

    assert cache.stats.canonical_solves==nsolve
    assert np.isclose(b.overlap,np.conj(a.overlap))
    assert np.allclose(b.centroid,np.conj(a.centroid))
    assert np.allclose(b.covariance,a.covariance)


def test_incremental_add_and_subset_match_full_rebuild():
    basis=_basis()
    provider=AnalyticCI2DFrameProvider(
        nuclear_mass_au=8.0
    )

    cache2=GaussianPairCache(basis[:2])
    S2,H2,N2=build_cached_spinor_lvc_matrices(
        cache2,provider
    )
    cache3=cache2.expanded(basis[2])
    S3,H3,N3=expand_cached_spinor_lvc_matrices(
        S2,H2,N2,cache3,provider
    )

    full_cache=GaussianPairCache(basis)
    Sf,Hf,Nf=build_cached_spinor_lvc_matrices(
        full_cache,provider
    )

    assert np.allclose(S3,Sf,atol=3e-13)
    assert np.allclose(H3,Hf,atol=3e-13)
    assert np.allclose(N3,Nf,atol=3e-13)

    # Old-old pairs were inherited; only N_new=3 new canonical pairs are solved.
    assert cache3.stats.inherited_pairs==3
    assert cache3.stats.canonical_solves==3

    Ssub,Hsub,Nsub,csub=subset_cached_spinor_lvc_matrices(
        S3,H3,N3,cache3,[0,2]
    )
    ref_cache=GaussianPairCache([basis[0],basis[2]])
    Sr,Hr,Nr=build_cached_spinor_lvc_matrices(
        ref_cache,provider
    )

    assert np.allclose(Ssub,Sr,atol=3e-13)
    assert np.allclose(Hsub,Hr,atol=3e-13)
    assert np.allclose(Nsub,Nr,atol=3e-13)
    assert csub.stats.inherited_pairs==3


def test_factorization_count_reduction_is_explicit():
    n=10
    assert v14_factorization_equivalent_for_sh(n)==385
    assert v14_factorization_equivalent_for_time(n)==300
    assert v15_factorization_count_for_snapshot(n)==55
