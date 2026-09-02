from dataclasses import dataclass, asdict
from pathlib import Path
import json
import time
import numpy as np

from .benchmark_campaign import CIPassageConfig
from .dynamic_gauge_graph import AnalyticCI2DFrameProvider
from .v16_benchmark import (
    _initial_problem,
    _exact_target,
    _chain_basis,
)
from .born_huang_grid_v12 import build_born_huang_grid_2d
from .exact2d import run_exact_2d
from .residual_basis_v13 import normalized_grid_density
from .electronic_observables import (
    density_matrix_populations,
    density_matrix_purity,
)
from .coherence_metrics import (
    coherence_phase_error,
    density_trace_distance,
)
from .error_controlled_sparse_dynamics_v17 import (
    ErrorControlledSparseSettingsV17,
    run_error_controlled_sparse_lvc_gaussians,
)
from .edge_importance_v17 import (
    EdgeImportanceSettingsV17,
    ErrorControlledGaussianLocalityGraphV17,
)
from .sparse_pair_matrices_v16 import (
    sparse_reduced_density,
    build_sparse_spinor_lvc_matrices,
)
from .pair_cache_v15 import (
    GaussianPairCache,
    build_cached_spinor_lvc_matrices,
)
from .sparse_error_budget_v17 import (
    score_threshold_snapshot_sweep,
    local_score_budget_snapshot_sweep,
    summarize_snapshot_convergence,
)


@dataclass(frozen=True)
class V17AcceptanceThresholds:
    max_initial_density_error: float = 0.035
    max_projected_dynamics_density_error: float = 0.001
    max_target_density_error: float = 0.035
    max_target_population_error: float = 0.03
    max_coherence_phase_error: float = 0.0035
    max_norm_drift: float = 1e-4
    max_condition_number: float = 5e3

    min_online_relaxations: int = 1
    max_final_audit_S_error: float = 0.006
    max_final_audit_H_error: float = 0.006
    max_final_audit_Snuc_error: float = 0.006
    max_unresolved_audits: int = 0
    max_local_omitted_score_l2: float = 0.0800001

    max_final_rho_difference_vs_v16: float = 5e-4

    max_finest_threshold_S_error: float = 0.001
    max_finest_threshold_H_error: float = 0.001

    min_n160_pair_reduction: float = 0.90
    max_active_edge_scaling_exponent: float = 1.15
    max_exact_pair_scaling_exponent: float = 1.20
    min_dense_pair_scaling_exponent: float = 1.90


def _complex_scalar(value):
    if isinstance(value,list) and len(value)==2:
        return complex(float(value[0]),float(value[1]))
    return complex(value)


def _complex_matrix(value):
    return np.asarray([
        [_complex_scalar(x) for x in row]
        for row in value
    ],dtype=complex)


def load_v16_context(repository_root):
    path=Path(repository_root)/"results"/"v016_sparse_locality_campaign.json"
    if not path.exists():
        return None
    data=json.loads(path.read_text(encoding="utf-8"))
    return {
        "reference":data["reference"],
        "final_density_matrix":
            _complex_matrix(
                data["reference"]["final_density_matrix"]
            ),
        "complexity":
            data["adaptive"]["complexity"],
        "final_sparse_matrix_audit":
            data["final_sparse_matrix_audit"],
    }


def run_edge_controlled_scaling_benchmark(
    sizes=(20,40,80,160),
    *,
    enter_score=0.03,
    exit_score=0.015,
    search_overlap_floor=1e-5,
    local_score_budget=0.08,
    dt=0.005,
):
    """Bounded-locality scaling benchmark for the actual v0.17 S/H/T scorer."""
    provider=AnalyticCI2DFrameProvider(
        nuclear_mass_au=20.0
    )
    rows=[]

    for n in sizes:
        basis=_chain_basis(n)
        graph=ErrorControlledGaussianLocalityGraphV17(
            provider,
            dt,
            EdgeImportanceSettingsV17(
                enter_score=enter_score,
                exit_score=exit_score,
                search_overlap_floor=
                    search_overlap_floor,
                overlap_weight=1.0,
                hamiltonian_weight=0.20,
                time_connection_weight=1.0,
                local_omitted_score_l2_budget=
                    local_score_budget,
            ),
        )

        t0=time.perf_counter()
        update=graph.update(basis)
        graph_seconds=time.perf_counter()-t0

        t0=time.perf_counter()
        sparse_mats=build_sparse_spinor_lvc_matrices(
            update,provider
        )
        sparse_seconds=time.perf_counter()-t0

        t0=time.perf_counter()
        dense_cache=GaussianPairCache(basis)
        build_cached_spinor_lvc_matrices(
            dense_cache,provider
        )
        dense_seconds=time.perf_counter()-t0

        dense_pairs=n*(n+1)//2
        actual=int(
            update.cache.stats.canonical_solves
        )

        rows.append({
            "n_basis":int(n),
            "active_edges":
                int(update.active_offdiagonal_edges),
            "total_offdiagonal_pairs":
                int(update.total_offdiagonal_pairs),
            "edge_fraction":
                float(update.edge_fraction),
            "sparsity_fraction":
                float(update.sparsity_fraction),
            "spatial_candidate_pairs":
                int(update.spatial_candidate_pairs),
            "pair_bound_screened_pairs":
                int(update.pair_bound_screened_pairs),
            "globally_screened_pairs":
                int(update.globally_screened_pairs),
            "exact_pair_checks":
                int(update.exact_pair_checks),
            "pair_factorizations":actual,
            "dense_canonical_pairs":
                int(dense_pairs),
            "pair_reduction_fraction":
                float(
                    1.0-actual/max(dense_pairs,1)
                ),
            "omitted_score_l2":
                float(update.omitted_candidate_score_l2),
            "budget_promoted_edges":
                int(update.budget_promoted_edges),
            "H_density":
                float(sparse_mats.H_density),
            "graph_seconds":
                float(graph_seconds),
            "sparse_matrix_seconds":
                float(sparse_seconds),
            "dense_matrix_seconds":
                float(dense_seconds),
            "assembly_speedup_vs_dense":
                float(
                    dense_seconds/max(
                        sparse_seconds,1e-30
                    )
                ),
        })

    return rows


def fit_v17_scaling_exponents(rows):
    N=np.asarray(
        [r["n_basis"] for r in rows],
        dtype=float,
    )

    def exponent(key):
        y=np.asarray(
            [r[key] for r in rows],
            dtype=float,
        )
        return float(
            np.polyfit(
                np.log(N),np.log(y),1
            )[0]
        )

    return {
        "active_edge_exponent":
            exponent("active_edges"),
        "spatial_candidate_exponent":
            exponent("spatial_candidate_pairs"),
        "exact_pair_check_exponent":
            exponent("exact_pair_checks"),
        "pair_factorization_exponent":
            exponent("pair_factorizations"),
        "dense_canonical_pair_exponent":
            exponent("dense_canonical_pairs"),
        "interpretation":(
            "Fitted exponents apply only to the deterministic bounded-locality "
            "Gaussian chain; worst-case dense configurations remain quadratic."
        ),
    }


def evaluate_v17_acceptance(
    reference,
    adaptive,
    final_audit,
    threshold_summary,
    scaling_rows,
    scaling_fit,
    v16_context=None,
    thresholds=None,
):
    t=thresholds or V17AcceptanceThresholds()

    unresolved=sum(
        1 for event in adaptive["events"]
        if event["kind"]=="sparse_audit_unresolved"
    )
    c=adaptive["complexity"]

    rho_diff=None
    if v16_context is not None:
        rho_diff=float(np.linalg.norm(
            reference["final_density_matrix"]
            -v16_context["final_density_matrix"],
            ord="fro",
        ))

    max_local_budget=max(
        [
            float(row["omitted_candidate_score_l2"])
            for row in adaptive["records"]
        ] or [0.0]
    )

    row160=next(
        row for row in scaling_rows
        if row["n_basis"]==160
    )

    checks={
        "initial_density_representation":
            reference["initial_density_error"]
            <=t.max_initial_density_error,
        "projected_dynamics":
            reference[
                "projected_dynamics_density_error"
            ]<=t.max_projected_dynamics_density_error,
        "target_density":
            reference["target_density_error"]
            <=t.max_target_density_error,
        "target_population":
            reference["target_population_error"]
            <=t.max_target_population_error,
        "coherence_phase":
            reference["coherence_phase_error"]
            is not None
            and reference["coherence_phase_error"]
            <=t.max_coherence_phase_error,
        "norm":
            reference["max_norm_drift"]
            <=t.max_norm_drift,
        "conditioning":
            reference["max_condition_number"]
            <=t.max_condition_number,

        "online_controller_relaxed":
            c["score_relaxations"]
            >=t.min_online_relaxations,
        "final_audit_S":
            final_audit[
                "relative_S_frobenius_error"
            ]<=t.max_final_audit_S_error,
        "final_audit_H":
            final_audit[
                "relative_H_frobenius_error"
            ]<=t.max_final_audit_H_error,
        "final_audit_Snuc":
            final_audit[
                "relative_Snuc_frobenius_error"
            ]<=t.max_final_audit_Snuc_error,
        "no_unresolved_audits":
            unresolved<=t.max_unresolved_audits,
        "local_importance_budget":
            max_local_budget
            <=t.max_local_omitted_score_l2,

        "rho_close_to_v16":
            rho_diff is None
            or rho_diff
            <=t.max_final_rho_difference_vs_v16,

        "threshold_S_converges":
            threshold_summary[
                "threshold_S_monotone"
            ],
        "threshold_H_converges":
            threshold_summary[
                "threshold_H_monotone"
            ],
        "budget_S_converges":
            threshold_summary[
                "budget_S_monotone"
            ],
        "budget_H_converges":
            threshold_summary[
                "budget_H_monotone"
            ],
        "finest_threshold_S":
            threshold_summary[
                "finest_threshold_S_error"
            ]<=t.max_finest_threshold_S_error,
        "finest_threshold_H":
            threshold_summary[
                "finest_threshold_H_error"
            ]<=t.max_finest_threshold_H_error,

        "n160_pair_reduction":
            row160["pair_reduction_fraction"]
            >=t.min_n160_pair_reduction,
        "active_edge_scaling":
            scaling_fit[
                "active_edge_exponent"
            ]<=t.max_active_edge_scaling_exponent,
        "exact_pair_scaling":
            scaling_fit[
                "exact_pair_check_exponent"
            ]<=t.max_exact_pair_scaling_exponent,
        "dense_pair_scaling":
            scaling_fit[
                "dense_canonical_pair_exponent"
            ]>=t.min_dense_pair_scaling_exponent,
    }

    return {
        "passed":bool(all(checks.values())),
        "checks":checks,
        "thresholds":asdict(t),
        "final_rho_difference_vs_v16":
            rho_diff,
        "maximum_recorded_local_omitted_score_l2":
            float(max_local_budget),
        "unresolved_audits":int(unresolved),
    }


def run_v017_release_benchmark(
    config=CIPassageConfig(),
    repository_root=None,
):
    provider=AnalyticCI2DFrameProvider(
        nuclear_mass_au=config.mass
    )
    initial_grid,psi_target,build=(
        _initial_problem(config,provider)
    )
    V,rho_target=_exact_target(
        config,provider,initial_grid,psi_target
    )

    exact_projected=run_exact_2d(
        build.projection.projected_wavefunction,
        initial_grid.dx,initial_grid.dx,V,
        mass=config.mass,
        dt=0.0025,
        steps=int(round(
            config.final_time/0.0025
        )),
        store_every=int(round(
            config.final_time/0.0025
        )),
    )
    rho_exact_projected=normalized_grid_density(
        exact_projected["psi"][-1],
        initial_grid.dx,
    )

    defect_grid=build_born_huang_grid_2d(
        grid_n=40,
        half_width=config.half_width,
        mass=config.mass,
        params=provider.params,
    )

    # Intentionally start too aggressive. The online dense audit must relax this
    # score threshold to the first setting satisfying the 0.6% S/H matrix budget.
    control=ErrorControlledSparseSettingsV17(
        defect_interval=10,
        enrich_relative_threshold=0.020,
        prune_relative_threshold=0.006,
        minimum_capture_fraction=0.003,
        minimum_local_utility=0.03,
        condition_penalty_weight=0.15,
        electronic_cost_weight=1.0,
        cost_horizon_steps=10,
        residual_shortlist=8,
        min_basis=8,
        max_basis=11,
        minimum_adaptation_separation_steps=10,
        minimum_prune_age_steps=20,
        prune_patience_checks=2,
        max_prune_fractional_loss=5e-7,
        max_replacement_prune_fractional_loss=5e-7,
        condition_limit=1e5,
        hard_condition_limit=5e6,
        candidate_position_shifts=(
            0.0,0.06,-0.06
        ),
        candidate_width_scales=(
            0.75,1.0,1.35
        ),
        edge_enter_score=0.060,
        edge_exit_score=0.030,
        search_overlap_floor=1e-5,
        edge_overlap_weight=1.0,
        edge_hamiltonian_weight=0.20,
        edge_time_connection_weight=1.0,
        local_omitted_score_l2_budget=0.08,
        sparsity_audit_interval=20,
        max_audit_S_error=0.006,
        max_audit_H_error=0.006,
        max_audit_Snuc_error=0.006,
        audit_relaxation_factor=0.5,
        max_audit_relaxations=4,
        check_initial_defect=False,
    )

    adaptive=run_error_controlled_sparse_lvc_gaussians(
        build.basis,
        C0=build.projection.coefficients,
        provider=provider,
        grid=defect_grid,
        dt=0.005,
        steps=int(round(
            config.final_time/0.005
        )),
        settings=control,
        store_every=10,
    )

    rho=sparse_reduced_density(
        adaptive["final_coefficients"],
        adaptive["final_sparse_matrices"].Snuc,
        normalize=True,
    )
    pop=density_matrix_populations(rho)
    target_pop=density_matrix_populations(
        rho_target
    )

    rho0_projected=normalized_grid_density(
        build.projection.projected_wavefunction,
        initial_grid.dx,
    )
    rho0_target=normalized_grid_density(
        psi_target,initial_grid.dx
    )

    reference={
        "initial_basis_size":10,
        "final_basis_size":
            len(adaptive["final_basis"]),
        "average_basis_size":
            adaptive["average_basis_size"],
        "projection_fidelity":
            float(build.projection.fidelity),
        "relative_residual":
            float(build.projection.relative_residual),
        "initial_density_error":float(
            np.linalg.norm(
                rho0_projected-rho0_target,
                ord="fro",
            )
        ),
        "projected_dynamics_density_error":
            float(np.linalg.norm(
                rho-rho_exact_projected,
                ord="fro",
            )),
        "target_density_error":float(
            np.linalg.norm(
                rho-rho_target,
                ord="fro",
            )
        ),
        "target_trace_distance":
            density_trace_distance(
                rho,rho_target
            ),
        "target_population_error":float(
            np.linalg.norm(
                pop-target_pop
            )
        ),
        "populations":pop,
        "target_populations":target_pop,
        "purity":
            density_matrix_purity(rho),
        "target_purity":
            density_matrix_purity(rho_target),
        "purity_error":float(abs(
            density_matrix_purity(rho)
            -density_matrix_purity(rho_target)
        )),
        "coherence":complex(rho[0,1]),
        "target_coherence":
            complex(rho_target[0,1]),
        "coherence_phase_error":
            coherence_phase_error(
                rho,rho_target
            ),
        "max_norm_drift":float(max(
            abs(row["norm"]-1.0)
            for row in adaptive["records"]
        )),
        "max_condition_number":float(max(
            row["condition_number"]
            for row in adaptive["records"]
        )),
        "final_density_matrix":rho,
    }

    final_audit=adaptive["audit_history"][-1]

    threshold_sweep=score_threshold_snapshot_sweep(
        adaptive["final_basis"],
        provider,
        dt=0.005,
        enter_scores=(
            0.12,0.08,0.06,0.04,
            0.03,0.02,0.01,
        ),
        search_overlap_floor=1e-5,
    )
    budget_sweep=local_score_budget_snapshot_sweep(
        adaptive["final_basis"],
        provider,
        dt=0.005,
        enter_score=0.06,
        budgets=(
            1e30,0.10,0.08,0.05,
            0.03,0.01,0.0,
        ),
        search_overlap_floor=1e-5,
    )
    sweep_summary=summarize_snapshot_convergence(
        threshold_sweep,budget_sweep
    )

    scaling=run_edge_controlled_scaling_benchmark()
    scaling_fit=fit_v17_scaling_exponents(
        scaling
    )

    v16=None
    if repository_root is not None:
        v16=load_v16_context(repository_root)

    acceptance=evaluate_v17_acceptance(
        reference,
        {
            "events":adaptive["events"],
            "records":adaptive["records"],
            "complexity":adaptive["complexity"],
        },
        final_audit,
        sweep_summary,
        scaling,
        scaling_fit,
        v16_context=v16,
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
        "reference":reference,
        "adaptive":{
            "events":adaptive["events"],
            "records":adaptive["records"],
            "defect_history":
                adaptive["defect_history"],
            "audit_history":
                adaptive["audit_history"],
            "cost_history":
                adaptive["cost_history"],
            "endpoint_graph":
                adaptive["endpoint_graph"],
            "midpoint_graph":
                adaptive["midpoint_graph"],
            "complexity":
                adaptive["complexity"],
            "control":asdict(control),
        },
        "final_audit":final_audit,
        "threshold_sweep":
            threshold_sweep,
        "budget_sweep":
            budget_sweep,
        "sweep_summary":
            sweep_summary,
        "scaling":
            scaling,
        "scaling_fit":
            scaling_fit,
        "v16_context":v16,
        "acceptance":acceptance,
    }
