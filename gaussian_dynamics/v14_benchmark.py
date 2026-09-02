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
    exact_reduced_electronic_density_diabatic,
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
from .adaptive_defect_dynamics_v14 import (
    AdaptiveDefectSettings,
    run_time_adaptive_defect_lvc_gaussians,
)
from .spinor_complete_lvc_v12 import (
    spinor_complete_reduced_density,
    coefficients_matrix,
)
from .residual_pruning_v14 import (
    prune_low_loss_gaussian_pair,
)
from .fast_lvc_matrices_v14 import (
    pair_evaluation_reduction,
)
from .complexity_v14 import asymptotic_complexity


@dataclass(frozen=True)
class V14AcceptanceThresholds:
    max_initial_density_error: float = 0.035
    max_projected_dynamics_density_error: float = 0.003
    max_target_density_error: float = 0.035
    max_target_population_error: float = 0.03
    max_coherence_phase_error: float = 0.0035
    max_norm_drift: float = 1e-4
    max_condition_number: float = 5e3
    min_enrichment_events: int = 1
    min_pair_evaluation_reduction: float = 0.40
    max_pruning_stress_loss: float = 1e-10


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
        node=("v14_seed",0),
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


def _pruning_stress_test(build,provider):
    """Add one nearly redundant zero-amplitude Gaussian and verify low-loss pruning."""
    basis=[
        DynamicGraphTBF(
            uid=int(b.uid),
            state=int(b.state),
            q=np.asarray(b.q,float).copy(),
            p=np.asarray(b.p,float).copy(),
            A=np.asarray(b.A,float).copy(),
            node=b.node,
        )
        for b in build.basis[:4]
    ]

    parent=basis[0]
    redundant=DynamicGraphTBF(
        uid=999999,
        state=int(parent.state),
        q=parent.q.copy(),
        p=parent.p.copy(),
        A=1.01*parent.A,
        node=("v14_pruning_stress",999999),
    )
    basis.append(redundant)

    from .fast_lvc_matrices_v14 import (
        build_spinor_complete_lvc_matrices_symmetric,
    )
    _,_,Snuc=build_spinor_complete_lvc_matrices_symmetric(
        basis,provider
    )

    Cfull=coefficients_matrix(
        build.projection.coefficients,
        len(build.basis),
    )
    Cmat=np.vstack([
        Cfull[:4],
        np.zeros((1,2),dtype=complex),
    ])

    result=prune_low_loss_gaussian_pair(
        Cmat,
        Snuc,
        uids=[b.uid for b in basis],
        max_fractional_loss=1e-8,
        protected_uids=[b.uid for b in basis[:-1]],
        require_condition_improvement=True,
    )
    if result is None:
        raise RuntimeError("v0.14 pruning stress test failed to find redundant pair.")

    return {
        "removed_uid":int(result.removed_uid),
        "fractional_projection_loss":
            float(result.fractional_projection_loss),
        "absolute_projection_loss":
            float(result.absolute_projection_loss),
        "condition_before":
            float(result.condition_before),
        "condition_after":
            float(result.condition_after),
        "condition_improvement_factor":
            float(result.condition_before/max(result.condition_after,1e-30)),
    }


def load_v13_context(repository_root):
    path=Path(repository_root)/"results"/"v013_residual_driven_campaign.json"
    if not path.exists():
        return None
    data=json.loads(path.read_text(encoding="utf-8"))
    r=data["reference"]
    return {
        "basis_size":r["basis_size"],
        "projection_fidelity":r["projection_fidelity"],
        "relative_residual":r["relative_residual"],
        "initial_density_error":r["initial_density_error"],
        "projected_dynamics_density_error":
            r["projected_dynamics_density_error"],
        "target_density_error":r["target_density_error"],
        "target_population_error":r["target_population_error"],
        "coherence_phase_error":r["coherence_phase_error"],
        "max_norm_drift":r["max_norm_drift"],
        "max_condition_number":r["max_condition_number"],
        "acceptance":data["acceptance"],
    }


def evaluate_v14_acceptance(
    reference,
    adaptive,
    pruning_stress,
    thresholds=None,
):
    t=thresholds or V14AcceptanceThresholds()

    enrichments=[
        e for e in adaptive["events"]
        if e["kind"]=="defect_enrichment"
    ]
    defect_reductions=all(
        e["relative_defect_after"]
        <e["relative_defect_before"]
        for e in enrichments
    )

    complexity=adaptive["complexity"]
    pair_reduction=1.0-(
        complexity["pair_matrix_evaluations"]
        /max(complexity["ordered_pair_equivalent"],1)
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
        "adaptive_enrichment":
            len(enrichments)>=t.min_enrichment_events,
        "enrichment_reduces_defect":
            bool(defect_reductions),
        "hermitian_pair_reduction":
            pair_reduction>=t.min_pair_evaluation_reduction,
        "low_loss_pruning":
            pruning_stress["fractional_projection_loss"]
            <=t.max_pruning_stress_loss,
        "pruning_improves_condition":
            pruning_stress["condition_after"]
            <pruning_stress["condition_before"],
    }

    return {
        "passed":bool(all(checks.values())),
        "checks":checks,
        "thresholds":asdict(t),
        "pair_evaluation_reduction":float(pair_reduction),
    }


def run_v014_release_benchmark(
    config=CIPassageConfig(),
    repository_root=None,
):
    provider=AnalyticCI2DFrameProvider(
        nuclear_mass_au=config.mass
    )
    initial_grid,psi_target,build=_initial_problem(
        config,provider
    )
    V,rho_target=_exact_target(
        config,provider,initial_grid,psi_target
    )

    # Exact propagation of the same projected state used by the adaptive Gaussian run.
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

    control=AdaptiveDefectSettings(
        defect_interval=10,
        enrich_relative_threshold=0.020,
        prune_relative_threshold=0.006,
        minimum_capture_fraction=0.003,
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

    adaptive=run_time_adaptive_defect_lvc_gaussians(
        build.basis,
        C0=build.projection.coefficients,
        provider=provider,
        grid=defect_grid,
        dt=0.005,
        steps=int(round(config.final_time/0.005)),
        settings=control,
        store_every=10,
    )
    rho=spinor_complete_reduced_density(
        adaptive["final_coefficients"],
        adaptive["final_basis"],
        normalize=True,
    )

    rho0_projected=normalized_grid_density(
        build.projection.projected_wavefunction,
        initial_grid.dx,
    )
    rho0_target=normalized_grid_density(
        psi_target,initial_grid.dx
    )

    target_pop=density_matrix_populations(
        rho_target
    )
    pop=density_matrix_populations(rho)

    reference={
        "initial_basis_size":10,
        "final_basis_size":len(adaptive["final_basis"]),
        "average_basis_size":adaptive["average_basis_size"],
        "projection_fidelity":
            float(build.projection.fidelity),
        "relative_residual":
            float(build.projection.relative_residual),
        "initial_density_error":float(np.linalg.norm(
            rho0_projected-rho0_target,ord="fro"
        )),
        "projected_dynamics_density_error":
            float(np.linalg.norm(
                rho-rho_exact_projected,ord="fro"
            )),
        "target_density_error":float(np.linalg.norm(
            rho-rho_target,ord="fro"
        )),
        "target_trace_distance":
            density_trace_distance(rho,rho_target),
        "target_population_error":float(np.linalg.norm(
            pop-target_pop
        )),
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
        "target_coherence":complex(rho_target[0,1]),
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

    pruning_stress=_pruning_stress_test(
        build,provider
    )
    acceptance=evaluate_v14_acceptance(
        reference,adaptive,pruning_stress
    )

    context=None
    if repository_root is not None:
        context=load_v13_context(repository_root)

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
            "defect_history":adaptive["defect_history"],
            "average_basis_size":adaptive["average_basis_size"],
            "complexity":adaptive["complexity"],
            "control":asdict(control),
        },
        "pruning_stress":pruning_stress,
        "complexity_model":
            asymptotic_complexity().__dict__.copy(),
        "half_build_pair_reduction_at_11":
            pair_evaluation_reduction(11),
        "acceptance":acceptance,
        "v13_context":context,
    }
