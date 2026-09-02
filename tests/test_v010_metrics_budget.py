import numpy as np

from gaussian_dynamics.error_budget import estimate_population_error_budget
from gaussian_dynamics.benchmark_acceptance import (
    BenchmarkThresholds,
    evaluate_managed_benchmark,
)


def fake_run():
    return {
        "records":[
            {
                "norm":1.0,
                "state_populations":np.array([0.2,0.8]),
                "condition_number":3.0,
                "basis_size":2,
                "spa1_relative_correction":0.01,
            },
            {
                "norm":1.0+1e-11,
                "state_populations":np.array([0.21,0.79]),
                "condition_number":4.0,
                "basis_size":2,
                "spa1_relative_correction":0.02,
            },
        ],
        "events":[
            {"kind":"spawn"},
            {"kind":"prune","projection_loss":1e-10},
        ],
    }


def test_acceptance_contract_passes_well_behaved_run():
    result=evaluate_managed_benchmark(
        fake_run(),
        reference_populations=np.array([0.2101,0.7899]),
        thresholds=BenchmarkThresholds(
            max_population_l2_vs_reference=1e-2
        ),
    )
    assert result.passed
    assert all(result.checks.values())


def test_error_budget_identifies_largest_sensitivity():
    budget=estimate_population_error_budget(
        exact_reference=[0.2,0.8],
        exact_next_coarser=[0.21,0.79],
        managed_reference_settings=[0.25,0.75],
        managed_next_coarser_dt=[0.26,0.74],
        spa0=[0.25,0.75],
        spa1=[0.20,0.80],
        spawn_threshold_low=[0.25,0.75],
        spawn_threshold_high=[0.251,0.749],
        basis_small=[0.25,0.75],
        basis_large=[0.249,0.751],
    )

    assert budget.spa_truncation_proxy > budget.managed_timestep_proxy
    assert budget.dominant_proxy=="spa_truncation_proxy"
