import time
import numpy as np

from .benchmark_campaign import CIPassageConfig
from .dynamic_gauge_graph import AnalyticCI2DFrameProvider
from .v16_benchmark import _initial_problem
from .born_huang_grid_v12 import build_born_huang_grid_2d
from .v18_benchmark import (
    ConvergenceCoordinatesV18,
    release_settings_v18,
    _exact_runs,
)
from .convergence_complete_dynamics_v18 import (
    run_convergence_complete_lvc_gaussians,
)
from .convergence_campaign_v18 import (
    evaluate_convergence_run_v18,
    compare_snapshot_trajectory_v18,
)
from .wavefunction_metrics_v18 import (
    gaussian_wavefunction_on_grid,
)


def _summarize_worker_run(
    coordinates,
    out,
    comparison_grid,
    exact_projected_final,
    exact_target_final,
    wall_seconds,
):
    metrics=evaluate_convergence_run_v18(
        out,
        comparison_grid,
        exact_projected_final,
        exact_target_final,
    )
    return {
        "coordinates":coordinates.as_dict(),
        "wall_seconds":float(wall_seconds),
        "basis_size":metrics["basis_size"],
        "average_basis_size":
            metrics["average_basis_size"],
        "projected_reduced_density_error":
            metrics["reduced_density_projected"][
                "density_frobenius_error"
            ],
        "target_reduced_density_error":
            metrics["reduced_density_target"][
                "density_frobenius_error"
            ],
        "target_population_error":
            metrics["reduced_density_target"][
                "population_l2_error"
            ],
        "target_coherence_phase_error":
            metrics["reduced_density_target"][
                "coherence_phase_error"
            ],
        "projected_wavefunction_fidelity":
            metrics["wavefunction_projected"][
                "fidelity"
            ],
        "projected_wavefunction_l2_error":
            metrics["wavefunction_projected"][
                "phase_aligned_l2_error"
            ],
        "projected_nuclear_density_l2_error":
            metrics["wavefunction_projected"][
                "nuclear_density_l2_error"
            ],
        "projected_nuclear_density_tv":
            metrics["wavefunction_projected"][
                "nuclear_density_total_variation"
            ],
        "projected_mean_error":
            metrics["wavefunction_projected"][
                "mean_error_l2"
            ],
        "projected_covariance_error":
            metrics["wavefunction_projected"][
                "covariance_error_frobenius"
            ],
        "target_wavefunction_fidelity":
            metrics["wavefunction_target"][
                "fidelity"
            ],
        "maximum_norm_drift":
            metrics["maximum_norm_drift"],
        "maximum_condition_number":
            metrics["maximum_condition_number"],
        "final_density_matrix":
            metrics["final_density_matrix"],
        "sentinel_audits":
            out["sentinel_audit_history"],
        "sampled_audits":
            out["sampled_audit_history"],
        "events":out["events"],
        "complexity":out["complexity"],
        "resolved_control_steps":
            out["settings"][
                "resolved_control_steps"
            ],
    }


def run_coordinate_worker_v18(
    coordinates,
    *,
    final_time=0.60,
    trajectory=False,
    trajectory_store_interval=0.10,
    include_final_wavefunction=False,
):
    """Self-contained convergence point intended for a fresh process.

    Gaussian propagation is executed before the exact-grid FFT reference. This avoids
    process-history effects from long FFT campaigns contaminating sparse-solver timing.
    """
    if isinstance(coordinates,dict):
        coordinates=ConvergenceCoordinatesV18(
            **coordinates
        )

    config=CIPassageConfig(
        final_time=float(final_time)
    )
    provider=AnalyticCI2DFrameProvider(
        nuclear_mass_au=config.mass
    )
    comparison_grid,psi_target,initial_build=(
        _initial_problem(config,provider)
    )
    defect_grid=build_born_huang_grid_2d(
        grid_n=40,
        half_width=config.half_width,
        mass=config.mass,
        params=provider.params,
    )

    settings=release_settings_v18(
        coordinates
    )
    steps=int(round(
        float(config.final_time)
        /float(coordinates.dt)
    ))
    store_every=max(
        int(round(
            float(trajectory_store_interval)
            /float(coordinates.dt)
        )),
        1,
    )

    t0=time.perf_counter()
    out=run_convergence_complete_lvc_gaussians(
        initial_build.basis,
        C0=initial_build.projection.coefficients,
        provider=provider,
        grid=defect_grid,
        dt=float(coordinates.dt),
        steps=steps,
        settings=settings,
        store_every=store_every,
        return_snapshots=bool(trajectory),
    )
    wall=time.perf_counter()-t0

    exact_target,exact_projected=_exact_runs(
        config,
        provider,
        comparison_grid,
        psi_target,
        initial_build.projection.projected_wavefunction,
        exact_dt=0.0025,
        store_interval_time=
            float(trajectory_store_interval),
    )

    row=_summarize_worker_run(
        coordinates,
        out,
        comparison_grid,
        exact_projected["psi"][-1],
        exact_target["psi"][-1],
        wall,
    )

    result={
        "result":row,
        "initial_projection":{
            "fidelity":
                float(initial_build.projection.fidelity),
            "relative_residual":
                float(initial_build.projection.relative_residual),
        },
    }

    if trajectory:
        result["trajectory"]=(
            compare_snapshot_trajectory_v18(
                out["snapshots"],
                exact_projected,
                comparison_grid,
            )
        )

    if include_final_wavefunction:
        result["final_wavefunction"]=(
            gaussian_wavefunction_on_grid(
                out["final_coefficients"],
                out["final_basis"],
                comparison_grid.points,
            )
        )
        result["comparison_grid_area"]=float(
            comparison_grid.dx*comparison_grid.dx
        )

    return result
