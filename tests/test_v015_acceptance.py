from gaussian_dynamics.v15_benchmark import (
    V15AcceptanceThresholds,
    evaluate_v15_acceptance,
)


def _reference():
    return {
        "initial_density_error":0.03,
        "projected_dynamics_density_error":1e-4,
        "target_density_error":0.033,
        "target_population_error":0.028,
        "coherence_phase_error":0.002,
        "max_norm_drift":1e-6,
        "max_condition_number":1500.0,
        "projection_fidelity":0.88,
        "purity":0.66,
    }


def _adaptive():
    return {
        "events":[{
            "kind":"cost_aware_defect_enrichment",
            "relative_defect_before":0.03,
            "relative_defect_after":0.025,
            "cost_aware_utility":0.2,
            "new_pair_factorizations_during_expansion":0,
        }],
        "complexity":{
            "factorization_reduction_fraction":0.85,
            "cache_hit_fraction":0.64,
        },
    }


def test_v15_acceptance_checks_cost_cache_and_physics_regression():
    v14={
        "reference":{
            "projection_fidelity":0.88,
            "initial_density_error":0.03,
            "projected_dynamics_density_error":1e-4,
            "target_density_error":0.033,
            "target_population_error":0.028,
            "purity":0.66,
            "coherence_phase_error":0.002,
            "max_norm_drift":1e-6,
            "max_condition_number":1500.0,
        }
    }

    out=evaluate_v15_acceptance(
        _reference(),_adaptive(),v14,
        V15AcceptanceThresholds(),
    )
    assert out["passed"]


def test_v15_acceptance_rejects_cache_regression():
    adaptive=_adaptive()
    adaptive["complexity"]["factorization_reduction_fraction"]=0.2

    out=evaluate_v15_acceptance(
        _reference(),adaptive,None
    )
    assert not out["passed"]
    assert not out["checks"]["pair_factorization_reduction"]
