from dataclasses import dataclass, asdict
import numpy as np


def _distance(a, b):
    return float(np.linalg.norm(np.asarray(a, float) - np.asarray(b, float)))


@dataclass(frozen=True)
class ErrorBudget:
    total_vs_exact: float
    exact_discretization_proxy: float
    managed_timestep_proxy: float
    spa_truncation_proxy: float
    spawn_threshold_proxy: float
    basis_size_proxy: float
    dominant_proxy: str

    def to_dict(self):
        return asdict(self)


def estimate_population_error_budget(
    exact_reference,
    exact_next_coarser,
    managed_reference_settings,
    managed_next_coarser_dt,
    spa0,
    spa1,
    spawn_threshold_low,
    spawn_threshold_high,
    basis_small,
    basis_large,
):
    """Construct a non-additive sensitivity/error budget from population vectors.

    These components are *proxies*, not statistically independent additive errors.
    They answer which controlled numerical/method change moves the result most.
    """
    values = {
        "exact_discretization_proxy": _distance(exact_reference, exact_next_coarser),
        "managed_timestep_proxy": _distance(
            managed_reference_settings, managed_next_coarser_dt
        ),
        "spa_truncation_proxy": _distance(spa0, spa1),
        "spawn_threshold_proxy": _distance(
            spawn_threshold_low, spawn_threshold_high
        ),
        "basis_size_proxy": _distance(basis_small, basis_large),
    }

    dominant = max(values, key=values.get)

    return ErrorBudget(
        total_vs_exact=_distance(managed_reference_settings, exact_reference),
        dominant_proxy=dominant,
        **values,
    )
