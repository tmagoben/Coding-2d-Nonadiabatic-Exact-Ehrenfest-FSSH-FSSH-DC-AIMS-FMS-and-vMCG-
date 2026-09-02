from dataclasses import dataclass, asdict
import numpy as np

from .initial_conditions import sample_gaussian_wigner
from .benchmark_campaign import (
    CIPassageConfig,
    run_exact_passage,
    run_managed_passage,
)


@dataclass(frozen=True)
class EnsembleStatistics:
    mean: np.ndarray
    std: np.ndarray
    sem: np.ndarray
    nsamples: int

    def to_dict(self):
        return {
            "mean": np.asarray(self.mean).tolist(),
            "std": np.asarray(self.std).tolist(),
            "sem": np.asarray(self.sem).tolist(),
            "nsamples": int(self.nsamples),
        }


def ensemble_statistics(values):
    x = np.asarray(values, dtype=float)
    if x.ndim < 2 or len(x) == 0:
        raise ValueError("values must have shape (nsample, ...).")

    n = len(x)
    mean = np.mean(x, axis=0)

    if n == 1:
        std = np.zeros_like(mean)
    else:
        std = np.std(x, axis=0, ddof=1)

    sem = std / np.sqrt(n)

    return EnsembleStatistics(mean=mean, std=std, sem=sem, nsamples=n)


def run_ci_initial_condition_ensemble(
    base_config=CIPassageConfig(),
    nsamples=4,
    seed=12345,
    exact_grid_n=48,
    exact_dt=0.005,
    managed_dt=0.005,
    spa_order=0,
    spawn_action_threshold=2e-4,
    max_basis=4,
):
    """Compare exact and managed dynamics over deterministic Gaussian Wigner samples."""
    samples = sample_gaussian_wigner(
        base_config.q_array(),
        base_config.p_array(),
        base_config.A_matrix(),
        nsamples=nsamples,
        seed=seed,
    )

    exact_populations = []
    managed_populations = []
    exact_norms = []
    managed_norms = []

    for q0, p0 in zip(samples.q, samples.p):
        config = CIPassageConfig(
            q0=tuple(float(x) for x in q0),
            p0=tuple(float(x) for x in p0),
            A_diag=base_config.A_diag,
            state=base_config.state,
            mass=base_config.mass,
            final_time=base_config.final_time,
            half_width=base_config.half_width,
        )

        exact = run_exact_passage(
            config,
            grid_n=exact_grid_n,
            dt=exact_dt,
        )
        managed = run_managed_passage(
            config,
            dt=managed_dt,
            spa_order=spa_order,
            spawn_action_threshold=spawn_action_threshold,
            max_basis=max_basis,
            store_every=max(1, int(round(config.final_time/managed_dt))),
        )

        exact_populations.append(exact["final_populations_adiabatic"])
        managed_populations.append(managed["records"][-1]["state_populations"])
        exact_norms.append(exact["norm"])
        managed_norms.append(managed["records"][-1]["norm"])

    exact_populations = np.asarray(exact_populations)
    managed_populations = np.asarray(managed_populations)
    differences = managed_populations - exact_populations

    return {
        "seed": int(seed),
        "samples_q": samples.q,
        "samples_p": samples.p,
        "exact_populations": exact_populations,
        "managed_populations": managed_populations,
        "population_differences": differences,
        "exact_statistics": ensemble_statistics(exact_populations),
        "managed_statistics": ensemble_statistics(managed_populations),
        "difference_statistics": ensemble_statistics(differences),
        "exact_norms": np.asarray(exact_norms),
        "managed_norms": np.asarray(managed_norms),
    }
