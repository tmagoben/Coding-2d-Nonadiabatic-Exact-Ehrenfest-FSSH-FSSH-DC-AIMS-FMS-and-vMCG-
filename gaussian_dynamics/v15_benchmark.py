from dataclasses import dataclass, asdict
from pathlib import Path
import json
import numpy as np

from .benchmark_campaign import CIPassageConfig
from .dynamic_gauge_graph import AnalyticCI2DFrameProvider
from .dynamic_graph_aims import DynamicGraphTBF
from .born_huang_grid_v12 import build_born_huang_grid_2d
from .exact_benchmark import localized_adiabatic_packet_2d
from .exact2d import run_exact_2d
from .ci2d import diabatic_potential_2d
from .electronic_observables import (
    density_matrix_populations,
    density_matrix_purity,
)
from .coherence_metrics import (
    coherence_phase_error,
    density_trace_distance,
)
from .residual_basis_v13 import (
    cartesian_offsets_2d,
    generate_gaussian_dictionary,
    prepare_gaussian_dictionary,
    build_residual_greedy_basis_prepared,
    normalized_grid_density,
)
from .adaptive_defect_dynamics_v15 import (
    AdaptiveDefectSettingsV15,
    run_time_adaptive_cost_aware_lvc_gaussians,
    reduced_density_from_snuc,
)
from .complexity_v15 import asymptotic_complexity_v15


@dataclass(frozen=True)
class V15AcceptanceThresholds:
    max_initial_density_error: float = 0.035
    max_projected_dynamics_density_error: float = 0.003
    max_target_density_error: float = 0.035
    max_target_population_error: float = 0.03
    max_coherence_phase_error: float = 0.0035
    max_norm_drift: float = 1e-4
    max_condition_number: float = 5e3

    min_enrichment_events: int = 1
    min_cost_aware_utility: float = 0.15
    min_factorization_reduction: float = 0.84
    min_cache_hit_fraction: float = 0.60
    max_incremental_expansion_pair_factorizations: int = 0

    max_v14_reference_metric_difference: float = 1e-9


def _initial_problem(config,provider):
    grid=build_born_huang_grid_2d(
        grid_n=64,
        half_width=config.half_width,
        mass=config.mass,
        params=provider.params,
    )
    psi=localized_adiabatic_packet_2d(
        grid.points,
        config.q_array(),
        config.p_array(),
        config.A_matrix(),
        state=config.state,
    )
    seed=DynamicGraphTBF(
        uid=0,
        state=config.state,
        q=config.q_array(),
        p=config.p_array(),
        A=config.A_matrix(),
        node=("v15_seed",0),
    )

    candidates=generate_gaussian_dictionary(
        config.q_array(),
        config.p_array(),
        config.A_matrix(),
        config.state,
        cartesian_offsets_2d(
            radius=1.0,
            spacing=0.2,
        ),
        width_scales=(
            1.0,1.5,2.0,3.0,4.0,6.0
        ),
    )
    prepared=prepare_gaussian_dictionary(
        candidates,grid.points,grid.dx
    )
    build=build_residual_greedy_basis_prepared(
        psi,
        grid.points,
        grid.dx,
        provider,
        [seed],
        prepared,
        max_basis=10,
        top_k_density_screen=30,
        density_screen=True,
        condition_limit=1e5,
    )
    return grid,psi,build


def _exact_target(config,provider,grid,psi):
    V=diabatic_potential_2d(
        grid.X,grid.Y,provider.params
    )
    exact=run_exact_2d(
        psi,
        grid.dx,grid.dx,V,
        mass=config.mass,
        dt=0.0025,
        steps=int(round(config.final_time/0.0025)),
        store_every=int(round(config.final_time/0.0025)),
    )
    rho=normalized_grid_density(
        exact["psi"][-1],grid.dx
    )
    return V,rho


def load_v14_context(repository_root):
    path=Path(repository_root)/"results"/"v014_time_adaptive_defect_campaign.json"
    if not path.exists():
        return None

    data=json.loads(
        path.read_text(encoding="utf-8")
    )
    r=data["reference"]
    c=data["adaptive"]["complexity"]
    e=[
        x for x in data["adaptive"]["events"]
        if x["kind"]=="defect_enrichment"
    ]
    return {
        "reference":{
            "initial_basis_size":r["initial_basis_size"],
            "final_basis_size":r["final_basis_size"],
            "average_basis_size":r["average_basis_size"],
            "projection_fidelity":r["projection_fidelity"],
            "initial_density_error":r["initial_density_error"],
            "projected_dynamics_density_error":
                r["projected_dynamics_density_error"],
            "target_density_error":r["target_density_error"],
            "target_population_error":r["target_population_error"],
            "purity":r["purity"],
            "coherence_phase_error":
                r["coherence_phase_error"],
            "max_norm_drift":r["max_norm_drift"],
            "max_condition_number":
                r["max_condition_number"],
        },
        "complexity":{
            "total_seconds":c["total_seconds"],
            "matrix_build_seconds":
                c["matrix_build_seconds"],
            "time_matrix_seconds":
                c["time_matrix_seconds"],
            "candidate_ranking_seconds":
                c["candidate_ranking_seconds"],
            "matrix_build_calls":
                c["matrix_build_calls"],
        },
        "event":None if not e else e[0],
    }


def _reference_metric_difference(current,v14):
    if v14 is None:
        return None

    keys=[
        "projection_fidelity",
        "initial_density_error",
        "projected_dynamics_density_error",
        "target_density_error",
        "target_population_error",
        "purity",
        "coherence_phase_error",
        "max_norm_drift",
        "max_condition_number",
    ]
    diffs={
        key:abs(
            float(current[key])
            -float(v14["reference"][key])
        )
        for key in keys
    }
    return {
        "per_metric":diffs,
        "maximum":float(max(diffs.values())),
    }


def evaluate_v15_acceptance(
    reference,
    adaptive,
    v14_context=None,
    thresholds=None,
):
    t=thresholds or V15AcceptanceThresholds()

    events=[
        e for e in adaptive["events"]
        if e["kind"]=="cost_aware_defect_enrichment"
    ]
    reductions=all(
        e["relative_defect_after"]
        <e["relative_defect_before"]
        for e in events
    )
    utilities=[
        e["cost_aware_utility"]
        for e in events
    ]
    max_expand_pair=max(
        [
            e["new_pair_factorizations_during_expansion"]
            for e in events
        ] or [0]
    )

    complexity=adaptive["complexity"]
    comparison=_reference_metric_difference(
        reference,v14_context
    )

    checks={
        "initial_density_representation":
            reference["initial_density_error"]
            <=t.max_initial_density_error,
        "projected_dynamics":
            reference["projected_dynamics_density_error"]
            <=t.max_projected_dynamics_density_error,
        "target_density":
            reference["target_density_error"]
            <=t.max_target_density_error,
        "target_population":
            reference["target_population_error"]
            <=t.max_target_population_error,
        "coherence_phase":
            reference["coherence_phase_error"] is not None
            and reference["coherence_phase_error"]
            <=t.max_coherence_phase_error,
        "norm":
            reference["max_norm_drift"]
            <=t.max_norm_drift,
        "conditioning":
            reference["max_condition_number"]
            <=t.max_condition_number,
        "cost_aware_enrichment":
            len(events)>=t.min_enrichment_events,
        "enrichment_reduces_defect":
            bool(reductions),
        "cost_utility_gate":
            bool(utilities)
            and min(utilities)>=t.min_cost_aware_utility,
        "pair_factorization_reduction":
            complexity["factorization_reduction_fraction"]
            >=t.min_factorization_reduction,
        "cache_reuse":
            complexity["cache_hit_fraction"]
            >=t.min_cache_hit_fraction,
        "incremental_expansion_reuses_candidate_pairs":
            max_expand_pair
            <=t.max_incremental_expansion_pair_factorizations,
        "v14_physics_regression":
            comparison is None
            or comparison["maximum"]
            <=t.max_v14_reference_metric_difference,
    }

    return {
        "passed":bool(all(checks.values())),
        "checks":checks,
        "thresholds":asdict(t),
        "v14_reference_difference":comparison,
    }


def run_v015_release_benchmark(
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
        steps=int(round(config.final_time/0.0025)),
        store_every=int(round(config.final_time/0.0025)),
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

    control=AdaptiveDefectSettingsV15(
        defect_interval=10,
        enrich_relative_threshold=0.020,
        prune_relative_threshold=0.006,
        minimum_capture_fraction=0.003,
        minimum_cost_aware_utility=0.15,
        condition_penalty_weight=0.15,
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
        candidate_position_shifts=(0.0,0.06,-0.06),
        candidate_width_scales=(0.75,1.0,1.35),
        check_initial_defect=False,
    )

    adaptive=run_time_adaptive_cost_aware_lvc_gaussians(
        build.basis,
        C0=build.projection.coefficients,
        provider=provider,
        grid=defect_grid,
        dt=0.005,
        steps=int(round(config.final_time/0.005)),
        settings=control,
        store_every=10,
    )

    rho=reduced_density_from_snuc(
        adaptive["final_coefficients"],
        adaptive["final_nuclear_overlap"],
        normalize=True,
    )

    rho0_projected=normalized_grid_density(
        build.projection.projected_wavefunction,
        initial_grid.dx,
    )
    rho0_target=normalized_grid_density(
        psi_target,initial_grid.dx
    )

    pop=density_matrix_populations(rho)
    target_pop=density_matrix_populations(
        rho_target
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
        "projected_dynamics_density_error":float(
            np.linalg.norm(
                rho-rho_exact_projected,
                ord="fro",
            )
        ),
        "target_density_error":float(
            np.linalg.norm(
                rho-rho_target,
                ord="fro",
            )
        ),
        "target_trace_distance":
            density_trace_distance(rho,rho_target),
        "target_population_error":float(
            np.linalg.norm(pop-target_pop)
        ),
        "populations":pop,
        "target_populations":target_pop,
        "purity":density_matrix_purity(rho),
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
            coherence_phase_error(rho,rho_target),
        "max_norm_drift":float(max(
            abs(r["norm"]-1.0)
            for r in adaptive["records"]
        )),
        "max_condition_number":float(max(
            r["condition_number_nuclear"]
            for r in adaptive["records"]
        )),
    }

    v14=None
    if repository_root is not None:
        v14=load_v14_context(repository_root)

    acceptance=evaluate_v15_acceptance(
        reference,adaptive,v14
    )

    timing_comparison=None
    if v14 is not None:
        old=float(
            v14["complexity"]["total_seconds"]
        )
        new=float(
            adaptive["complexity"]["total_seconds"]
        )
        timing_comparison={
            "v14_saved_adaptive_seconds":old,
            "v15_adaptive_seconds":new,
            "saved_benchmark_speedup":
                float(old/max(new,1e-30)),
            "saved_benchmark_runtime_reduction":
                float(1.0-new/max(old,1e-30)),
            "note":(
                "Wall times are environment dependent and are reported diagnostically; "
                "they are not release acceptance criteria."
            ),
        }

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
            "defect_history":
                adaptive["defect_history"],
            "cost_history":
                adaptive["cost_history"],
            "average_basis_size":
                adaptive["average_basis_size"],
            "complexity":
                adaptive["complexity"],
            "control":asdict(control),
        },
        "v14_context":v14,
        "timing_comparison":timing_comparison,
        "complexity_model":
            asymptotic_complexity_v15().__dict__.copy(),
        "acceptance":acceptance,
    }
