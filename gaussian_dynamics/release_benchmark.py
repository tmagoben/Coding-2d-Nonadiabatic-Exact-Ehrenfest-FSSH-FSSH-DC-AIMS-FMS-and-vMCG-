import numpy as np

from .benchmark_campaign import (
    CIPassageConfig,
    run_exact_grid_timestep_surface,
    select_finest_exact_reference,
    run_exact_passage,
    run_managed_passage,
)
from .electronic_observables import (
    exact_reduced_electronic_density_diabatic,
    reduced_electronic_density_analytic_ci_diabatic,
    density_matrix_populations,
    density_matrix_purity,
    density_matrix_linear_entropy,
)
from .benchmark_metrics import summarize_managed_run
from .benchmark_acceptance import (
    BenchmarkThresholds,
    evaluate_managed_benchmark,
)
from .error_budget import estimate_population_error_budget


def _normalized_exact_diabatic_density(out):
    rho = exact_reduced_electronic_density_diabatic(
        out["psi_final"], out["dx"], out["dx"]
    )
    return rho/np.trace(rho)


def _managed_observables(run):
    rho = reduced_electronic_density_analytic_ci_diabatic(
        run["final_coefficients"],
        run["final_basis"],
        normalize=True,
    )
    return {
        "rho": rho,
        "populations": density_matrix_populations(rho),
        "purity": density_matrix_purity(rho),
        "linear_entropy": density_matrix_linear_entropy(rho),
    }


def run_compact_v010_release_benchmark(
    config=CIPassageConfig(),
):
    """Deterministic compact campaign used as a release regression report."""
    exact_surface = run_exact_grid_timestep_surface(
        config,
        grid_values=(32, 48, 64),
        dt_values=(0.010, 0.005, 0.0025),
    )
    finest_row = select_finest_exact_reference(exact_surface)

    # Re-run the candidate reference because the surface intentionally stores only
    # compact final observables rather than the full final wavefunction.
    exact_ref = run_exact_passage(config, grid_n=64, dt=0.0025)
    exact_next = run_exact_passage(config, grid_n=48, dt=0.005)

    rho_exact = _normalized_exact_diabatic_density(exact_ref)
    rho_exact_next = _normalized_exact_diabatic_density(exact_next)

    p_exact = density_matrix_populations(rho_exact)
    p_exact_next = density_matrix_populations(rho_exact_next)

    common = dict(
        config=config,
        spawn_action_threshold=2e-4,
        overlap_block=0.9999,
        minimum_spawn_separation_steps=5,
        store_every=20,
    )

    runs = {
        "reference": run_managed_passage(
            dt=0.005, spa_order=0, max_basis=4, **common
        ),
        "coarse_dt": run_managed_passage(
            dt=0.010, spa_order=0, max_basis=4, **common
        ),
        "spa1": run_managed_passage(
            dt=0.005, spa_order=1, max_basis=4, **common
        ),
        "spawn_refined": run_managed_passage(
            dt=0.005,
            spa_order=0,
            max_basis=4,
            **{**common, "spawn_action_threshold":1e-4},
        ),
        "basis_small": run_managed_passage(
            dt=0.005, spa_order=0, max_basis=2, **common
        ),
        "basis_large": run_managed_passage(
            dt=0.005, spa_order=0, max_basis=6, **common
        ),
    }

    obs = {name:_managed_observables(run) for name,run in runs.items()}

    budget = estimate_population_error_budget(
        exact_reference=p_exact,
        exact_next_coarser=p_exact_next,
        managed_reference_settings=obs["reference"]["populations"],
        managed_next_coarser_dt=obs["coarse_dt"]["populations"],
        spa0=obs["reference"]["populations"],
        spa1=obs["spa1"]["populations"],
        spawn_threshold_low=obs["spawn_refined"]["populations"],
        spawn_threshold_high=obs["reference"]["populations"],
        basis_small=obs["basis_small"]["populations"],
        basis_large=obs["basis_large"]["populations"],
    )

    acceptance = evaluate_managed_benchmark(
        runs["reference"],
        reference_populations=p_exact,
        observed_populations=obs["reference"]["populations"],
        thresholds=BenchmarkThresholds(
            max_norm_error=1e-6,
            max_population_sum_error=1e-8,
            max_condition_number=1e8,
            max_total_pruning_loss=1e-6,
            max_population_l2_vs_reference=5e-2,
        ),
    )

    compact_runs = {}
    for name,run in runs.items():
        metrics = summarize_managed_run(run)
        compact_runs[name] = {
            "metrics": metrics.to_dict(),
            "diabatic_populations": obs[name]["populations"],
            "purity": obs[name]["purity"],
            "linear_entropy": obs[name]["linear_entropy"],
            "events": run["events"],
        }

    return {
        "config": {
            "q0": list(config.q0),
            "p0": list(config.p0),
            "A_diag": list(config.A_diag),
            "state": config.state,
            "mass": config.mass,
            "final_time": config.final_time,
            "half_width": config.half_width,
        },
        "exact_surface": exact_surface,
        "candidate_exact_reference": finest_row,
        "exact_diabatic_populations": p_exact,
        "exact_purity": density_matrix_purity(rho_exact),
        "exact_linear_entropy": density_matrix_linear_entropy(rho_exact),
        "managed_runs": compact_runs,
        "error_budget": budget.to_dict(),
        "acceptance": acceptance.to_dict(),
    }
