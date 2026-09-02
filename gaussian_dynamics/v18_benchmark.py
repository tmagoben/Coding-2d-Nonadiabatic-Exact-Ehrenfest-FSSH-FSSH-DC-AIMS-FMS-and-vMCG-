from dataclasses import dataclass, asdict, replace
from pathlib import Path
import json
import numpy as np

from .benchmark_campaign import CIPassageConfig
from .dynamic_gauge_graph import AnalyticCI2DFrameProvider
from .v16_benchmark import _initial_problem, _chain_basis
from .born_huang_grid_v12 import build_born_huang_grid_2d
from .ci2d import diabatic_potential_2d
from .exact2d import run_exact_2d
from .residual_basis_v13 import normalized_grid_density
from .convergence_complete_dynamics_v18 import (
    ConvergenceCompleteSettingsV18,
    run_convergence_complete_lvc_gaussians,
)
from .convergence_campaign_v18 import (
    ConvergenceCoordinatesV18,
    evaluate_convergence_run_v18,
    compare_snapshot_trajectory_v18,
    successive_self_convergence_order,
)
from .wavefunction_metrics_v18 import (
    gaussian_wavefunction_on_grid,
    compare_wavefunctions,
    phase_aligned_l2_error,
)
from .sparse_pair_matrices_v16 import (
    sparse_reduced_density,
)
from .edge_importance_v17 import (
    EdgeImportanceSettingsV17,
    ErrorControlledGaussianLocalityGraphV17,
)
from .sampled_sparse_audit_v18 import (
    sampled_omitted_edge_audit_v18,
)


@dataclass(frozen=True)
class V18AcceptanceThresholds:
    min_projected_fidelity: float = 0.980
    max_projected_phase_l2: float = 0.140
    max_projected_density_l2: float = 0.060
    max_projected_mean_error: float = 0.003
    max_projected_covariance_error: float = 0.015

    max_reduced_projected_error: float = 2e-4
    max_reduced_target_error: float = 0.035
    max_norm_drift: float = 1e-4
    max_condition_number: float = 1e4

    min_trajectory_fidelity: float = 0.975
    max_exact_projection_fidelity_drift: float = 1e-9

    max_sampled_audit_failures: int = 0
    max_sentinel_failures: int = 0
    required_dense_sentinels: int = 2
    min_dense_audit_pair_reduction_vs_v17: float = 0.65
    min_candidate_peak_memory_reduction: float = 0.95

    min_basis_ladder_improvement_fraction: float = 0.25
    max_dt_medium_fine_self_l2: float = 1e-3
    min_dt_self_convergence_order: float = 1.5
    min_edge_budget_improvement: float = 1e-4
    min_growth_trigger_improvement_fraction: float = 0.20


def _complex_scalar(value):
    if isinstance(value,list) and len(value)==2:
        return complex(float(value[0]),float(value[1]))
    return complex(value)


def _complex_matrix(value):
    return np.asarray([
        [_complex_scalar(x) for x in row]
        for row in value
    ],dtype=complex)


def load_v17_context(repository_root):
    path=Path(repository_root)/"results"/"v017_sparse_error_control_campaign.json"
    if not path.exists():
        return None
    data=json.loads(path.read_text(encoding="utf-8"))
    return {
        "reference":data["reference"],
        "complexity":data["adaptive"]["complexity"],
        "final_density_matrix":
            _complex_matrix(
                data["reference"]["final_density_matrix"]
            ),
    }


def _base_settings():
    return ConvergenceCompleteSettingsV18(
        defect_interval=10,
        defect_interval_time=0.05,
        enrich_relative_threshold=0.015,
        prune_relative_threshold=0.004,
        minimum_capture_fraction=0.002,
        minimum_local_utility=0.01,
        condition_penalty_weight=0.15,
        electronic_cost_weight=1.0,
        cost_horizon_steps=10,
        cost_horizon_time=0.05,
        residual_shortlist=8,

        min_basis=8,
        max_basis=13,
        minimum_adaptation_separation_steps=10,
        minimum_adaptation_separation_time=0.05,
        minimum_prune_age_steps=20,
        minimum_prune_age_time=0.10,
        prune_patience_checks=2,

        max_prune_fractional_loss=5e-7,
        max_replacement_prune_fractional_loss=5e-7,
        emergency_prune_fractional_loss=1e-4,

        condition_limit=1e5,
        hard_condition_limit=5e6,
        orthogonal_norm_floor=1e-8,

        candidate_position_shifts=(0.0,0.06,-0.06),
        candidate_width_scales=(0.75,1.0,1.35),
        candidate_momentum_directions=("nac","momentum"),
        include_same_surface_candidates=True,
        include_other_surface_candidates=True,
        candidate_overlap_block=0.999999,

        edge_enter_score=0.030,
        edge_exit_score=0.015,
        search_overlap_floor=5e-6,
        edge_overlap_weight=1.0,
        edge_hamiltonian_weight=0.20,
        edge_time_connection_weight=1.0,
        local_omitted_score_l2_budget=0.010,

        sampled_audit_interval=20,
        sampled_audit_interval_time=0.10,
        sampled_audit_priority_pairs=8,
        sampled_audit_random_pairs=8,
        sampled_audit_wider_search_factor=0.1,
        sampled_audit_seed=20260813,
        sampled_audit_violation_factor=1.0,
        sampled_audit_relaxation_factor=0.5,
        max_sampled_audit_relaxations=3,

        sentinel_max_S_error=0.006,
        sentinel_max_H_error=0.006,
        sentinel_max_Snuc_error=0.006,

        candidate_batch_size=16,
        check_initial_defect=False,
    )



def release_settings_v18(coordinates):
    """Return the canonical v0.18 controls for one convergence coordinate."""
    if isinstance(coordinates,dict):
        coordinates=ConvergenceCoordinatesV18(
            **coordinates
        )
    return replace(
        _base_settings(),
        max_basis=int(coordinates.max_basis),
        local_omitted_score_l2_budget=
            float(coordinates.local_score_budget),
        enrich_relative_threshold=
            float(coordinates.enrich_threshold),
    )


def _exact_runs(
    config,
    provider,
    comparison_grid,
    psi_target,
    projected_wavefunction,
    *,
    exact_dt=0.0025,
    store_interval_time=0.10,
):
    """Exact target/projected grid trajectories with matched storage cadence."""
    V=diabatic_potential_2d(
        comparison_grid.X,
        comparison_grid.Y,
        provider.params,
    )
    steps=int(round(
        float(config.final_time)/float(exact_dt)
    ))
    store_every=max(
        int(round(
            float(store_interval_time)
            /float(exact_dt)
        )),
        1,
    )
    target=run_exact_2d(
        psi_target,
        comparison_grid.dx,
        comparison_grid.dx,
        V,
        mass=config.mass,
        dt=float(exact_dt),
        steps=steps,
        store_every=store_every,
    )
    projected=run_exact_2d(
        projected_wavefunction,
        comparison_grid.dx,
        comparison_grid.dx,
        V,
        mass=config.mass,
        dt=float(exact_dt),
        steps=steps,
        store_every=store_every,
    )
    return target,projected


def sampled_audit_scaling_v18(
    sizes=(20,40,80,160),
    *,
    dt=0.005,
    priority_count=8,
    random_count=8,
):
    """Bounded-locality audit-cost diagnostic.

    The graph construction still performs exact local S/H/T scoring. This diagnostic
    isolates the *normal audit* cost: sampled omitted-edge checks remain bounded by the
    requested sample count rather than rebuilding every dense pair.
    """
    provider=AnalyticCI2DFrameProvider(
        nuclear_mass_au=20.0
    )
    rows=[]
    for n in sizes:
        basis=_chain_basis(int(n))
        settings=EdgeImportanceSettingsV17(
            enter_score=0.030,
            exit_score=0.015,
            search_overlap_floor=5e-6,
            overlap_weight=1.0,
            hamiltonian_weight=0.20,
            time_connection_weight=1.0,
            local_omitted_score_l2_budget=0.010,
        )
        graph=ErrorControlledGaussianLocalityGraphV17(
            provider,float(dt),settings
        )
        update=graph.update(basis)
        audit=sampled_omitted_edge_audit_v18(
            basis,provider,float(dt),
            update,graph.settings,
            step=20,
            priority_count=int(priority_count),
            random_count=int(random_count),
            seed=20260813,
        )
        dense_pairs=int(n)*(int(n)-1)//2
        rows.append({
            "n_basis":int(n),
            "active_edges":
                int(update.active_offdiagonal_edges),
            "omitted_pairs":
                int(audit.omitted_pair_count_estimate),
            "sampled_pairs":
                int(len(audit.sampled_pairs)),
            "all_offdiagonal_pairs":
                int(dense_pairs),
            "sample_fraction_of_all_pairs":
                float(
                    len(audit.sampled_pairs)
                    /max(dense_pairs,1)
                ),
            "maximum_sampled_score":
                float(audit.maximum_score),
            "passed":bool(audit.passed),
        })
    return rows


def assemble_v018_campaign_from_partials(partials):
    """Assemble fresh-process worker outputs into a compact coordinate table.

    This utility is intentionally simple: it accepts an iterable of worker result
    dictionaries returned by `run_coordinate_worker_v18` and sorts them by
    `(dt, max_basis, local_score_budget, enrich_threshold)`.
    """
    rows=[]
    trajectories=[]
    initial_projection=None
    for item in partials:
        if "result" not in item:
            raise ValueError("each partial must contain a 'result' entry.")
        rows.append(item["result"])
        if "trajectory" in item:
            trajectories.append({
                "coordinates":
                    item["result"]["coordinates"],
                "trajectory":
                    item["trajectory"],
            })
        if initial_projection is None:
            initial_projection=item.get(
                "initial_projection"
            )
    rows.sort(key=lambda row:(
        float(row["coordinates"]["dt"]),
        int(row["coordinates"]["max_basis"]),
        float(row["coordinates"]["local_score_budget"]),
        float(row["coordinates"]["enrich_threshold"]),
    ))
    return {
        "rows":rows,
        "trajectories":trajectories,
        "initial_projection":initial_projection,
    }

def _reference_problem(config,provider):
    comparison_grid,psi_target,build=_initial_problem(
        config,provider
    )
    V=diabatic_potential_2d(
        comparison_grid.X,
        comparison_grid.Y,
        provider.params,
    )

    exact_dt=0.0025
    store_dt=0.10
    store_every=int(round(store_dt/exact_dt))
    steps=int(round(config.final_time/exact_dt))

    exact_target=run_exact_2d(
        psi_target,
        comparison_grid.dx,
        comparison_grid.dx,
        V,
        mass=config.mass,
        dt=exact_dt,
        steps=steps,
        store_every=store_every,
    )
    exact_projected=run_exact_2d(
        build.projection.projected_wavefunction,
        comparison_grid.dx,
        comparison_grid.dx,
        V,
        mass=config.mass,
        dt=exact_dt,
        steps=steps,
        store_every=store_every,
    )

    area=float(
        comparison_grid.dx*comparison_grid.dx
    )
    overlap_rows=[]
    for t,psi_t,psi_p in zip(
        exact_target["time"],
        exact_target["psi"],
        exact_projected["psi"],
    ):
        metrics=compare_wavefunctions(
            psi_t,psi_p,
            comparison_grid.points,
            area,
        )
        overlap_rows.append({
            "time":float(t),
            "fidelity":float(metrics["fidelity"]),
            "phase_aligned_l2_error":
                float(metrics["phase_aligned_l2_error"]),
        })

    fidelity=np.asarray([
        row["fidelity"]
        for row in overlap_rows
    ])
    exact_overlap_summary={
        "rows":overlap_rows,
        "initial_fidelity":float(fidelity[0]),
        "final_fidelity":float(fidelity[-1]),
        "maximum_fidelity_drift":float(
            np.max(np.abs(fidelity-fidelity[0]))
        ),
        "interpretation":(
            "Exact target and exact projected states evolve under the same unitary "
            "Hamiltonian, so their overlap should remain constant. This separates "
            "initial representation error from Gaussian propagation error."
        ),
    }

    return (
        comparison_grid,
        psi_target,
        build,
        exact_target,
        exact_projected,
        exact_overlap_summary,
    )


def _run_key(dt,settings):
    return (
        round(float(dt),12),
        int(settings.max_basis),
        round(
            float(settings.local_omitted_score_l2_budget),
            12,
        ),
        round(
            float(settings.enrich_relative_threshold),
            12,
        ),
    )


def _jsonable_metrics(metrics):
    """Drop arrays/complex phase factors that are not needed in axis tables."""
    return {
        "fidelity":float(metrics["fidelity"]),
        "phase_aligned_l2_error":
            float(metrics["phase_aligned_l2_error"]),
        "nuclear_density_l2_error":
            float(metrics["nuclear_density_l2_error"]),
        "nuclear_density_total_variation":
            float(metrics["nuclear_density_total_variation"]),
        "mean_error_l2":
            float(metrics["mean_error_l2"]),
        "covariance_error_frobenius":
            float(metrics["covariance_error_frobenius"]),
    }


def evaluate_v18_acceptance(
    canonical,
    trajectory,
    exact_overlap,
    basis_axis,
    dt_axis,
    dt_self,
    edge_axis,
    growth_axis,
    v17_context,
    thresholds=None,
):
    t=thresholds or V18AcceptanceThresholds()
    proj=canonical["wavefunction_projected"]
    reduced_proj=canonical[
        "reduced_density_projected"
    ]
    reduced_target=canonical[
        "reduced_density_target"
    ]
    complexity=canonical["complexity"]

    sentinel_failures=sum(
        1 for row in canonical["sentinel_audits"]
        if not row["passed"]
    )

    dense_pair_reduction=None
    rho_diff_v17=None
    if v17_context is not None:
        old_pairs=float(
            v17_context["complexity"][
                "audit_pair_factorizations"
            ]
        )
        new_pairs=float(
            complexity[
                "sentinel_pair_factorizations"
            ]
        )
        dense_pair_reduction=float(
            1.0-new_pairs/max(old_pairs,1.0)
        )
        rho_diff_v17=float(np.linalg.norm(
            canonical["final_density_matrix"]
            -v17_context["final_density_matrix"],
            ord="fro",
        ))

    basis_errors=[
        row["wavefunction_projected"][
            "phase_aligned_l2_error"
        ]
        for row in basis_axis
    ]
    basis_monotone=all(
        b<a
        for a,b in zip(
            basis_errors[:-1],
            basis_errors[1:],
        )
    )
    basis_improvement=float(
        (basis_errors[0]-basis_errors[-1])
        /max(abs(basis_errors[0]),1e-30)
    )

    dt_self_errors=[
        row["phase_aligned_l2_error"]
        for row in dt_self
    ]
    dt_monotone=all(
        b<a
        for a,b in zip(
            dt_self_errors[:-1],
            dt_self_errors[1:],
        )
    )

    edge_errors=[
        row["wavefunction_projected"][
            "phase_aligned_l2_error"
        ]
        for row in edge_axis
    ]
    edge_nonincreasing=all(
        b<=a+1e-12
        for a,b in zip(
            edge_errors[:-1],
            edge_errors[1:],
        )
    )
    edge_improvement=float(
        edge_errors[0]-edge_errors[-1]
    )

    growth_errors=[
        row["wavefunction_projected"][
            "phase_aligned_l2_error"
        ]
        for row in growth_axis
    ]
    growth_nonincreasing=all(
        b<=a+1e-12
        for a,b in zip(
            growth_errors[:-1],
            growth_errors[1:],
        )
    )
    growth_improvement=float(
        (growth_errors[0]-growth_errors[-1])
        /max(abs(growth_errors[0]),1e-30)
    )

    min_traj_fidelity=float(min(
        row["fidelity"] for row in trajectory
    ))

    checks={
        "projected_fidelity":
            proj["fidelity"]
            >=t.min_projected_fidelity,
        "projected_phase_l2":
            proj["phase_aligned_l2_error"]
            <=t.max_projected_phase_l2,
        "projected_density_l2":
            proj["nuclear_density_l2_error"]
            <=t.max_projected_density_l2,
        "projected_mean_error":
            proj["mean_error_l2"]
            <=t.max_projected_mean_error,
        "projected_covariance_error":
            proj["covariance_error_frobenius"]
            <=t.max_projected_covariance_error,
        "reduced_projected":
            reduced_proj["density_frobenius_error"]
            <=t.max_reduced_projected_error,
        "reduced_target":
            reduced_target["density_frobenius_error"]
            <=t.max_reduced_target_error,
        "norm":
            canonical["maximum_norm_drift"]
            <=t.max_norm_drift,
        "conditioning":
            canonical["maximum_condition_number"]
            <=t.max_condition_number,

        "trajectory_fidelity":
            min_traj_fidelity
            >=t.min_trajectory_fidelity,
        "exact_projection_overlap_conserved":
            exact_overlap["maximum_fidelity_drift"]
            <=t.max_exact_projection_fidelity_drift,

        "sampled_audits":
            complexity["sampled_audit_failures"]
            <=t.max_sampled_audit_failures,
        "sentinel_failures":
            sentinel_failures
            <=t.max_sentinel_failures,
        "two_dense_sentinels":
            complexity["sentinel_dense_audits"]
            ==t.required_dense_sentinels,
        "dense_audit_pair_reduction_vs_v17":
            dense_pair_reduction is None
            or dense_pair_reduction
            >=t.min_dense_audit_pair_reduction_vs_v17,
        "candidate_peak_memory":
            complexity[
                "candidate_peak_memory_reduction_fraction"
            ]>=t.min_candidate_peak_memory_reduction,

        "basis_ladder_monotone":
            basis_monotone,
        "basis_ladder_improves":
            basis_improvement
            >=t.min_basis_ladder_improvement_fraction,

        "dt_self_monotone":
            dt_monotone,
        "dt_medium_fine":
            dt_self_errors[-1]
            <=t.max_dt_medium_fine_self_l2,
        "dt_self_order":
            dt_axis["observed_self_order"] is not None
            and dt_axis["observed_self_order"]
            >=t.min_dt_self_convergence_order,

        "edge_budget_nonincreasing":
            edge_nonincreasing,
        "edge_budget_improves":
            edge_improvement
            >=t.min_edge_budget_improvement,

        "growth_trigger_nonincreasing":
            growth_nonincreasing,
        "growth_trigger_improves":
            growth_improvement
            >=t.min_growth_trigger_improvement_fraction,
    }

    return {
        "passed":bool(all(checks.values())),
        "checks":checks,
        "thresholds":asdict(t),
        "minimum_trajectory_fidelity":
            min_traj_fidelity,
        "dense_audit_pair_reduction_vs_v17":
            dense_pair_reduction,
        "final_rho_difference_vs_v17":
            rho_diff_v17,
        "basis_ladder_improvement_fraction":
            basis_improvement,
        "edge_budget_absolute_improvement":
            edge_improvement,
        "growth_trigger_improvement_fraction":
            growth_improvement,
    }


def run_v018_release_benchmark(
    config=CIPassageConfig(),
    repository_root=None,
):
    provider=AnalyticCI2DFrameProvider(
        nuclear_mass_au=config.mass
    )
    (
        comparison_grid,
        psi_target,
        build,
        exact_target,
        exact_projected,
        exact_overlap,
    )=_reference_problem(config,provider)

    defect_grid=build_born_huang_grid_2d(
        grid_n=40,
        half_width=config.half_width,
        mass=config.mass,
        params=provider.params,
    )

    base=_base_settings()
    cache={}
    metrics_cache={}
    wavefunction_cache={}

    def run_one(
        *,
        dt=0.005,
        max_basis=None,
        local_budget=None,
        enrich_threshold=None,
        snapshots=False,
    ):
        settings=base
        if max_basis is not None:
            settings=replace(
                settings,max_basis=int(max_basis)
            )
        if local_budget is not None:
            settings=replace(
                settings,
                local_omitted_score_l2_budget=
                    float(local_budget),
            )
        if enrich_threshold is not None:
            settings=replace(
                settings,
                enrich_relative_threshold=
                    float(enrich_threshold),
            )

        key=_run_key(dt,settings)
        if key not in cache:
            steps=int(round(
                config.final_time/float(dt)
            ))
            store_every=(
                max(
                    int(round(0.10/float(dt))),
                    1,
                )
                if snapshots
                else steps
            )
            cache[key]=run_convergence_complete_lvc_gaussians(
                build.basis,
                C0=build.projection.coefficients,
                provider=provider,
                grid=defect_grid,
                dt=float(dt),
                steps=steps,
                settings=settings,
                store_every=store_every,
                return_snapshots=bool(snapshots),
            )
        elif snapshots and not cache[key]["snapshots"]:
            # Re-run the canonical coordinate only when trajectory snapshots are
            # explicitly requested after a final-only cached call.
            steps=int(round(
                config.final_time/float(dt)
            ))
            cache[key]=run_convergence_complete_lvc_gaussians(
                build.basis,
                C0=build.projection.coefficients,
                provider=provider,
                grid=defect_grid,
                dt=float(dt),
                steps=steps,
                settings=settings,
                store_every=max(
                    int(round(0.10/float(dt))),
                    1,
                ),
                return_snapshots=True,
            )

        if key not in metrics_cache:
            metrics_cache[key]=evaluate_convergence_run_v18(
                cache[key],
                comparison_grid,
                exact_projected["psi"][-1],
                exact_target["psi"][-1],
            )
            wavefunction_cache[key]=gaussian_wavefunction_on_grid(
                cache[key]["final_coefficients"],
                cache[key]["final_basis"],
                comparison_grid.points,
            )
        return key,settings,cache[key],metrics_cache[key]

    # Canonical release coordinate.
    canonical_key,canonical_settings,canonical_run,canonical_metrics=run_one(
        dt=0.005,
        max_basis=13,
        local_budget=0.010,
        enrich_threshold=0.015,
        snapshots=True,
    )
    canonical_trajectory=compare_snapshot_trajectory_v18(
        canonical_run["snapshots"],
        exact_projected,
        comparison_grid,
    )

    # Basis-completeness ladder.
    basis_axis=[]
    for nmax in (10,11,12,13):
        key,settings,out,metrics=run_one(
            dt=0.005,
            max_basis=nmax,
            local_budget=0.010,
            enrich_threshold=0.015,
        )
        basis_axis.append({
            "max_basis":int(nmax),
            "final_basis_size":
                int(metrics["basis_size"]),
            "average_basis_size":
                float(metrics["average_basis_size"]),
            "wavefunction_projected":
                _jsonable_metrics(
                    metrics["wavefunction_projected"]
                ),
            "reduced_density_projected":
                metrics["reduced_density_projected"],
            "enrichment_steps":[
                int(e["step"])
                for e in out["events"]
                if e["kind"]
                =="sparse_cost_aware_enrichment"
            ],
        })

    # Timestep convergence with all adaptive controls expressed in physical time.
    dt_rows=[]
    for dt in (0.010,0.005,0.0025):
        key,settings,out,metrics=run_one(
            dt=dt,
            max_basis=13,
            local_budget=0.010,
            enrich_threshold=0.015,
        )
        dt_rows.append({
            "dt":float(dt),
            "final_basis_size":
                int(metrics["basis_size"]),
            "wavefunction_projected":
                _jsonable_metrics(
                    metrics["wavefunction_projected"]
                ),
            "reduced_density_projected":
                metrics["reduced_density_projected"],
            "resolved_control_steps":
                out["settings"][
                    "resolved_control_steps"
                ],
            "_key":key,
        })

    area=float(
        comparison_grid.dx*comparison_grid.dx
    )
    dt_self=[]
    for a,b in zip(dt_rows[:-1],dt_rows[1:]):
        psi_a=wavefunction_cache[a["_key"]]
        psi_b=wavefunction_cache[b["_key"]]
        err=phase_aligned_l2_error(
            psi_b,psi_a,area
        )
        dt_self.append({
            "dt_coarse":float(a["dt"]),
            "dt_fine":float(b["dt"]),
            "phase_aligned_l2_error":
                float(err),
        })
    self_order=successive_self_convergence_order(
        dt_self[0]["phase_aligned_l2_error"],
        dt_self[1]["phase_aligned_l2_error"],
        refinement_ratio=2.0,
    )
    for row in dt_rows:
        row.pop("_key",None)
    dt_axis={
        "rows":dt_rows,
        "successive_solution_differences":
            dt_self,
        "observed_self_order":self_order,
        "interpretation":(
            "Observed order is computed from successive Gaussian solution differences, "
            "not from error to the exact grid state. This isolates timestep "
            "discretization from the nonzero Gaussian basis/model error floor."
        ),
    }

    # Sparse edge-budget convergence.
    edge_axis=[]
    for budget in (0.030,0.010,0.0):
        key,settings,out,metrics=run_one(
            dt=0.005,
            max_basis=13,
            local_budget=budget,
            enrich_threshold=0.015,
        )
        edge_axis.append({
            "local_score_budget":
                float(budget),
            "final_basis_size":
                int(metrics["basis_size"]),
            "average_graph_sparsity":
                float(
                    metrics["complexity"][
                        "average_sparsity_fraction"
                    ]
                ),
            "wavefunction_projected":
                _jsonable_metrics(
                    metrics["wavefunction_projected"]
                ),
            "final_sentinel_H_error":
                float(
                    out["sentinel_audit_history"][-1][
                        "relative_H_frobenius_error"
                    ]
                ),
        })

    # Adaptive growth-trigger sensitivity.
    growth_axis=[]
    for threshold in (
        0.050,0.035,0.030,0.025,0.015
    ):
        key,settings,out,metrics=run_one(
            dt=0.005,
            max_basis=13,
            local_budget=0.010,
            enrich_threshold=threshold,
        )
        growth_axis.append({
            "enrich_relative_threshold":
                float(threshold),
            "final_basis_size":
                int(metrics["basis_size"]),
            "enrichment_steps":[
                int(e["step"])
                for e in out["events"]
                if e["kind"]
                =="sparse_cost_aware_enrichment"
            ],
            "wavefunction_projected":
                _jsonable_metrics(
                    metrics["wavefunction_projected"]
                ),
        })

    v17=None
    if repository_root is not None:
        v17=load_v17_context(repository_root)

    canonical={
        **canonical_metrics,
        "wavefunction_projected":
            _jsonable_metrics(
                canonical_metrics[
                    "wavefunction_projected"
                ]
            ),
        "wavefunction_target":
            _jsonable_metrics(
                canonical_metrics[
                    "wavefunction_target"
                ]
            ),
        "sentinel_audits":
            canonical_run[
                "sentinel_audit_history"
            ],
        "sampled_audits":
            canonical_run[
                "sampled_audit_history"
            ],
        "events":
            canonical_run["events"],
        "settings":{
            "dt":0.005,
            "control":asdict(
                canonical_settings
            ),
            "resolved_control_steps":
                canonical_run["settings"][
                    "resolved_control_steps"
                ],
        },
    }

    acceptance=evaluate_v18_acceptance(
        canonical,
        canonical_trajectory,
        exact_overlap,
        basis_axis,
        dt_axis,
        dt_self,
        edge_axis,
        growth_axis,
        v17,
    )

    return {
        "config":{
            "q0":list(config.q0),
            "p0":list(config.p0),
            "A_diag":list(config.A_diag),
            "state":config.state,
            "mass":config.mass,
            "final_time":config.final_time,
            "half_width":config.half_width,
        },
        "initial_representation":{
            "basis_size":10,
            "projection_fidelity":
                float(build.projection.fidelity),
            "relative_residual":
                float(build.projection.relative_residual),
        },
        "exact_projected_target_overlap":
            exact_overlap,
        "canonical":canonical,
        "trajectory_projected_reference":
            canonical_trajectory,
        "basis_axis":basis_axis,
        "dt_axis":dt_axis,
        "edge_budget_axis":edge_axis,
        "growth_trigger_axis":growth_axis,
        "v17_context":v17,
        "run_cache_size":int(len(cache)),
        "acceptance":acceptance,
    }
