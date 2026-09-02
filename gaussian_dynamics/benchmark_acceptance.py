from dataclasses import dataclass, asdict
import numpy as np

from .benchmark_metrics import summarize_managed_run, population_sum_error


@dataclass(frozen=True)
class BenchmarkThresholds:
    max_norm_error: float = 1e-8
    max_population_sum_error: float = 1e-8
    max_condition_number: float = 1e10
    max_total_pruning_loss: float = 1e-6
    max_population_l2_vs_reference: float = 5e-3


@dataclass(frozen=True)
class BenchmarkAcceptance:
    passed: bool
    checks: dict
    metrics: dict

    def to_dict(self):
        return {
            "passed": self.passed,
            "checks": dict(self.checks),
            "metrics": dict(self.metrics),
        }


def evaluate_managed_benchmark(
    run,
    reference_populations=None,
    observed_populations=None,
    thresholds=None,
):
    """Evaluate propagation checks and, optionally, one explicitly supplied observable.

    `observed_populations` should be used when the rigorous population observable is
    not the legacy TBF-state-label proxy stored in the run records.
    """
    thresholds = thresholds or BenchmarkThresholds()
    m = summarize_managed_run(run, None)

    populations = (
        np.asarray(observed_populations, dtype=float)
        if observed_populations is not None
        else m.final_populations
    )

    checks = {
        "norm": m.max_norm_error <= thresholds.max_norm_error,
        "population_sum": population_sum_error(populations)
        <= thresholds.max_population_sum_error,
        "conditioning": m.max_condition_number <= thresholds.max_condition_number,
        "pruning_loss": m.total_pruning_loss <= thresholds.max_total_pruning_loss,
    }

    exact_error = None
    if reference_populations is not None:
        reference = np.asarray(reference_populations, dtype=float)
        exact_error = float(np.linalg.norm(populations-reference))
        checks["exact_reference_population"] = (
            exact_error <= thresholds.max_population_l2_vs_reference
        )

    metrics = m.to_dict()
    metrics["observed_populations"] = populations.tolist()
    metrics["observed_population_l2_vs_reference"] = exact_error

    return BenchmarkAcceptance(
        passed=bool(all(checks.values())),
        checks=checks,
        metrics=metrics,
    )
