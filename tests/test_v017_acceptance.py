import numpy as np

from gaussian_dynamics.v17_benchmark import (
    V17AcceptanceThresholds,
    evaluate_v17_acceptance,
)


def test_v17_acceptance_requires_relaxation_audit_and_convergence():
    reference={
        "initial_density_error":0.03,
        "projected_dynamics_density_error":2e-4,
        "target_density_error":0.033,
        "target_population_error":0.028,
        "coherence_phase_error":0.002,
        "max_norm_drift":1e-6,
        "max_condition_number":1200.0,
        "final_density_matrix":
            np.array([[0.2,0.1],[0.1,0.8]],complex),
    }
    adaptive={
        "events":[],
        "records":[
            {"omitted_candidate_score_l2":0.05}
        ],
        "complexity":{
            "score_relaxations":1,
        },
    }
    final_audit={
        "relative_S_frobenius_error":0.004,
        "relative_H_frobenius_error":0.004,
        "relative_Snuc_frobenius_error":0.004,
    }
    summary={
        "threshold_S_monotone":True,
        "threshold_H_monotone":True,
        "budget_S_monotone":True,
        "budget_H_monotone":True,
        "finest_threshold_S_error":2e-4,
        "finest_threshold_H_error":2e-4,
    }
    scaling=[
        {
            "n_basis":160,
            "pair_reduction_fraction":0.94,
        }
    ]
    fit={
        "active_edge_exponent":1.05,
        "exact_pair_check_exponent":1.08,
        "dense_canonical_pair_exponent":1.98,
    }
    v16={
        "final_density_matrix":
            np.array([[0.2,0.1],[0.1,0.8]],complex)
    }

    out=evaluate_v17_acceptance(
        reference,adaptive,final_audit,
        summary,scaling,fit,v16,
        V17AcceptanceThresholds(),
    )
    assert out["passed"]


def test_v17_acceptance_rejects_unresolved_audit():
    reference={
        "initial_density_error":0.03,
        "projected_dynamics_density_error":2e-4,
        "target_density_error":0.033,
        "target_population_error":0.028,
        "coherence_phase_error":0.002,
        "max_norm_drift":1e-6,
        "max_condition_number":1200.0,
        "final_density_matrix":np.eye(2)/2,
    }
    adaptive={
        "events":[
            {"kind":"sparse_audit_unresolved"}
        ],
        "records":[
            {"omitted_candidate_score_l2":0.05}
        ],
        "complexity":{
            "score_relaxations":1,
        },
    }
    audit={
        "relative_S_frobenius_error":0.004,
        "relative_H_frobenius_error":0.004,
        "relative_Snuc_frobenius_error":0.004,
    }
    summary={
        "threshold_S_monotone":True,
        "threshold_H_monotone":True,
        "budget_S_monotone":True,
        "budget_H_monotone":True,
        "finest_threshold_S_error":2e-4,
        "finest_threshold_H_error":2e-4,
    }
    scaling=[{
        "n_basis":160,
        "pair_reduction_fraction":0.94,
    }]
    fit={
        "active_edge_exponent":1.05,
        "exact_pair_check_exponent":1.08,
        "dense_canonical_pair_exponent":1.98,
    }

    out=evaluate_v17_acceptance(
        reference,adaptive,audit,
        summary,scaling,fit,None
    )
    assert not out["passed"]
    assert not out["checks"]["no_unresolved_audits"]
