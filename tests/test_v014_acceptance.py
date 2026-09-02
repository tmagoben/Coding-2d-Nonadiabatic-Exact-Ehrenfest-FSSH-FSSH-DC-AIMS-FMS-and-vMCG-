import numpy as np

from gaussian_dynamics.v14_benchmark import (
    V14AcceptanceThresholds,
    evaluate_v14_acceptance,
)


def test_v14_acceptance_checks_adaptation_and_complexity():
    reference={
        "initial_density_error":0.02,
        "projected_dynamics_density_error":0.001,
        "target_density_error":0.025,
        "target_population_error":0.02,
        "coherence_phase_error":0.002,
        "max_norm_drift":1e-6,
        "max_condition_number":1000.0,
    }
    adaptive={
        "events":[{
            "kind":"defect_enrichment",
            "relative_defect_before":0.04,
            "relative_defect_after":0.03,
        }],
        "complexity":{
            "pair_matrix_evaluations":55,
            "ordered_pair_equivalent":100,
        },
    }
    pruning={
        "fractional_projection_loss":1e-12,
        "condition_before":1e6,
        "condition_after":1e3,
    }

    result=evaluate_v14_acceptance(
        reference,adaptive,pruning,
        V14AcceptanceThresholds(),
    )
    assert result["passed"]
    assert np.isclose(
        result["pair_evaluation_reduction"],0.45
    )


def test_v14_acceptance_rejects_nonreducing_enrichment():
    reference={
        "initial_density_error":0.02,
        "projected_dynamics_density_error":0.001,
        "target_density_error":0.025,
        "target_population_error":0.02,
        "coherence_phase_error":0.002,
        "max_norm_drift":1e-6,
        "max_condition_number":1000.0,
    }
    adaptive={
        "events":[{
            "kind":"defect_enrichment",
            "relative_defect_before":0.03,
            "relative_defect_after":0.04,
        }],
        "complexity":{
            "pair_matrix_evaluations":55,
            "ordered_pair_equivalent":100,
        },
    }
    pruning={
        "fractional_projection_loss":0.0,
        "condition_before":1e6,
        "condition_after":1e3,
    }

    result=evaluate_v14_acceptance(
        reference,adaptive,pruning
    )
    assert not result["passed"]
    assert not result["checks"]["enrichment_reduces_defect"]
