import numpy as np

from gaussian_dynamics.analytic_molecular_backend_v19 import (
    AnalyticMolecularLVCBackendV19,
    default_diatomic_two_mode_map_v19,
)
from gaussian_dynamics.indexed_molecular_provider_v20 import (
    IndexedTrackedMolecularDirectProviderV20,
)
from gaussian_dynamics.sparse_molecular_matrices_v20 import (
    SparseMolecularTBFV20,
    MolecularSparseSettingsV20,
    SparseMolecularEdgeGraphV20,
    build_sparse_molecular_matrices_v20,
    build_dense_molecular_reference_v20,
)


def _basis(n=5,spacing=0.8):
    A=1.4*np.eye(2)
    x0=-0.5*spacing*(n-1)
    return [
        SparseMolecularTBFV20(
            uid=i,
            state=i%2,
            q=np.array([x0+i*spacing,0.35]),
            p=np.array([0.1*(-1)**i,0.0]),
            A=A,
        )
        for i in range(n)
    ]


def _provider():
    gmap=default_diatomic_two_mode_map_v19()
    return IndexedTrackedMolecularDirectProviderV20(
        AnalyticMolecularLVCBackendV19(gmap),
        gmap,
        rebuild_batch=4,
    )


def test_sparse_builder_equals_dense_when_all_edges_retained():
    provider=_provider()
    basis=_basis(5,0.7)
    settings=MolecularSparseSettingsV20(
        enter_score=1e-14,
        exit_score=1e-14,
        search_overlap_floor=1e-14,
        local_omitted_score_l2_budget=0.0,
        use_kdtree=False,
    )
    graph=SparseMolecularEdgeGraphV20(
        provider,0.01,settings
    )
    update=graph.update(basis)
    mats=build_sparse_molecular_matrices_v20(
        basis,update
    )
    dense=build_dense_molecular_reference_v20(
        basis,provider,0.01,settings
    )

    assert len(update.active_edges)==10
    assert np.allclose(
        mats.S.toarray(),dense["S"],atol=2e-12
    )
    assert np.allclose(
        mats.H.toarray(),dense["H"],atol=2e-12
    )
    assert np.allclose(
        mats.T_seed.toarray(),dense["T_seed"],atol=2e-12
    )


def test_sparse_candidate_centroids_are_subquadratic_on_local_chain():
    provider=_provider()
    basis=_basis(30,2.0)
    settings=MolecularSparseSettingsV20(
        enter_score=0.03,
        exit_score=0.015,
        search_overlap_floor=1e-4,
        local_omitted_score_l2_budget=0.01,
    )
    graph=SparseMolecularEdgeGraphV20(
        provider,0.01,settings
    )
    update=graph.update(basis)

    total=30*29//2
    assert update.exact_pair_checks<0.25*total
    assert len(update.active_edges)<0.25*total
