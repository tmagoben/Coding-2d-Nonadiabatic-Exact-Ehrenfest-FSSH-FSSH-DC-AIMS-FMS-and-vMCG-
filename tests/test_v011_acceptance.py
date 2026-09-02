from gaussian_dynamics.v11_benchmark import (
    V11AcceptanceThresholds,
    evaluate_v11_acceptance,
)


def test_v11_acceptance_distinguishes_population_from_density():
    metrics={
        "population_l2_error":0.02,
        "density_frobenius_error":0.15,
        "purity_error":0.03,
        "max_norm_drift":1e-3,
        "max_condition_number":1e4,
    }

    result=evaluate_v11_acceptance(
        metrics,
        V11AcceptanceThresholds(
            max_population_l2_error=0.05,
            max_density_frobenius_error=0.10,
            max_purity_error=0.05,
            max_norm_drift=0.01,
            max_condition_number=1e6,
        ),
    )

    assert result["checks"]["population"] is True
    assert result["checks"]["full_density"] is False
    assert result["passed"] is False
