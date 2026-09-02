from gaussian_dynamics.v16_benchmark import (
    V16AcceptanceThresholds,
    evaluate_v16_acceptance,
)


def test_v16_acceptance_requires_real_sparsity_and_dense_reference_agreement():
    reference={
        "initial_density_error":0.03,
        "projected_dynamics_density_error":2e-4,
        "target_density_error":0.033,
        "target_population_error":0.028,
        "coherence_phase_error":0.002,
        "max_norm_drift":1e-6,
        "max_condition_number":1500.0,
        "final_density_matrix":[[0.2,0.1],[0.1,0.8]],
    }
    adaptive={
        "events":[],
        "complexity":{
            "average_sparsity_fraction":0.10,
            "propagation_pair_factorizations":100,
        },
    }
    scaling=[
        {
            "n_basis":20,
            "pair_reduction_fraction":0.7,
            "edge_fraction":0.2,
        },
        {
            "n_basis":40,
            "pair_reduction_fraction":0.85,
            "edge_fraction":0.1,
        },
        {
            "n_basis":80,
            "pair_reduction_fraction":0.94,
            "edge_fraction":0.05,
        },
    ]
    demo={
        "cached_geometry":{
            "normalized_incremental_cost":0.4,
        },
        "new_geometry":{
            "normalized_incremental_cost":2.0,
        },
    }
    v15={
        "final_density_matrix":
            __import__("numpy").array([[0.2,0.1],[0.1,0.8]],complex),
        "complexity":{
            "propagation_pair_factorizations":200,
        },
    }

    out=evaluate_v16_acceptance(
        reference,adaptive,scaling,demo,v15,
        V16AcceptanceThresholds(),
    )
    assert out["passed"]


def test_v16_acceptance_rejects_graph_that_is_effectively_dense():
    import numpy as np

    reference={
        "initial_density_error":0.03,
        "projected_dynamics_density_error":2e-4,
        "target_density_error":0.033,
        "target_population_error":0.028,
        "coherence_phase_error":0.002,
        "max_norm_drift":1e-6,
        "max_condition_number":1500.0,
        "final_density_matrix":np.eye(2)/2,
    }
    adaptive={
        "events":[],
        "complexity":{
            "average_sparsity_fraction":0.001,
            "propagation_pair_factorizations":100,
        },
    }
    scaling=[
        {"n_basis":80,"pair_reduction_fraction":0.94,"edge_fraction":0.05}
    ]
    demo={
        "cached_geometry":{"normalized_incremental_cost":0.4},
        "new_geometry":{"normalized_incremental_cost":2.0},
    }

    out=evaluate_v16_acceptance(
        reference,adaptive,scaling,demo,None
    )
    assert not out["passed"]
    assert not out["checks"]["graph_is_actually_sparse"]
