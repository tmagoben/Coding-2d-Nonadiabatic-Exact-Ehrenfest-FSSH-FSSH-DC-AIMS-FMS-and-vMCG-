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
)
from gaussian_dynamics.sparse_molecular_dynamics_v20 import (
    SparseMolecularDynamicsSettingsV20,
    run_sparse_molecular_dynamics_v20,
    run_dense_molecular_reference_dynamics_v20,
)


def _provider():
    gmap=default_diatomic_two_mode_map_v19()
    return IndexedTrackedMolecularDirectProviderV20(
        AnalyticMolecularLVCBackendV19(gmap),
        gmap,
        rebuild_batch=4,
    )


def _basis():
    A=1.2*np.eye(2)
    return [
        SparseMolecularTBFV20(
            0,0,np.array([-0.8,0.35]),
            np.array([1.0,0.0]),A
        ),
        SparseMolecularTBFV20(
            1,1,np.array([0.0,0.38]),
            np.array([0.2,0.1]),A
        ),
        SparseMolecularTBFV20(
            2,0,np.array([0.8,0.34]),
            np.array([-0.8,0.0]),A
        ),
    ]


def test_sparse_all_edge_dynamics_matches_dense_reference():
    graph=MolecularSparseSettingsV20(
        enter_score=1e-14,
        exit_score=1e-14,
        search_overlap_floor=1e-14,
        local_omitted_score_l2_budget=0.0,
        use_kdtree=False,
    )
    settings=SparseMolecularDynamicsSettingsV20(
        graph=graph,
        sampled_audit_interval=100,
        dense_sentinel_S_limit=1e-11,
        dense_sentinel_H_limit=1e-11,
        dense_sentinel_T_limit=1e-11,
    )

    C0=np.array([
        1.0+0.0j,
        0.1+0.04j,
        -0.05j,
    ])
    sparse_out=run_sparse_molecular_dynamics_v20(
        _basis(),C0,_provider(),
        dt=0.002,steps=12,
        settings=settings,store_every=3,
    )
    dense_out=run_dense_molecular_reference_dynamics_v20(
        _basis(),C0,_provider(),
        dt=0.002,steps=12,store_every=3,
    )

    assert np.allclose(
        sparse_out["final_coefficients"],
        dense_out["final_coefficients"],
        atol=2e-10,
    )
    for a,b in zip(
        sparse_out["final_basis"],
        dense_out["final_basis"],
    ):
        assert np.allclose(a.q,b.q,atol=2e-12)
        assert np.allclose(a.p,b.p,atol=2e-12)

    assert sparse_out["sentinels"]["initial"]["passed"]
    assert sparse_out["sentinels"]["final"]["passed"]
    norms=np.asarray([
        row["norm"] for row in sparse_out["records"]
    ])
    assert np.max(np.abs(norms-1.0))<2e-8


def test_sampled_audit_relaxes_overaggressive_geometric_search():
    basis=_basis()[:2]
    graph=MolecularSparseSettingsV20(
        enter_score=0.03,
        exit_score=0.015,
        search_overlap_floor=0.90,
        local_omitted_score_l2_budget=0.01,
    )
    settings=SparseMolecularDynamicsSettingsV20(
        graph=graph,
        sampled_audit_interval=1,
        sampled_audit_priority_pairs=2,
        sampled_audit_random_pairs=0,
        sampled_audit_search_factor=0.1,
        sampled_audit_relaxation_factor=0.5,
        max_sampled_audit_relaxations=3,
        dense_sentinel_S_limit=1.0,
        dense_sentinel_H_limit=1.0,
        dense_sentinel_T_limit=1.0,
    )
    out=run_sparse_molecular_dynamics_v20(
        basis,
        np.array([1.0+0j,0.05+0j]),
        _provider(),
        dt=0.001,
        steps=1,
        settings=settings,
        store_every=1,
    )

    relax=[
        e for e in out["events"]
        if e["kind"]=="sampled_molecular_search_relaxation"
    ]
    assert len(relax)>=1
    assert relax[0]["new_search_overlap_floor"]<0.90
    assert out["sampled_audits"][-1]["passed"]
