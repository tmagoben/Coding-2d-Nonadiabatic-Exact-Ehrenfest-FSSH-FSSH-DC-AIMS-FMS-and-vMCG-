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
    coherence_magnitude,
    coherence_magnitude_error,
    coherence_phase_error,
    density_trace_distance,
    bloch_vector_error,
)
from .spinor_complete_dynamics_v12 import (
    run_spinor_complete_lvc_gaussians,
)
from .spinor_complete_lvc_v12 import (
    spinor_complete_reduced_density,
)
from .residual_basis_v13 import (
    cartesian_offsets_2d,
    generate_gaussian_dictionary,
    build_residual_greedy_basis,
    build_residual_greedy_basis_prepared,
    prepare_gaussian_dictionary,
    normalized_grid_density,
)
from .tdse_defect_v13 import (
    enrich_basis_from_tdse_defect,
)


@dataclass(frozen=True)
class V13AcceptanceThresholds:
    max_initial_density_error: float = 0.033
    max_projected_dynamics_density_error: float = 2e-4
    max_target_density_error: float = 0.033
    max_target_population_error: float = 0.03
    max_coherence_phase_error: float = 0.003
    max_norm_drift: float = 1e-4
    max_condition_number: float = 5e3
    min_defect_squared_reduction: float = 1e-8
    max_defect_prediction_relative_error: float = 5e-3


def _seed(config):
    return DynamicGraphTBF(
        uid=0,
        state=config.state,
        q=config.q_array(),
        p=config.p_array(),
        A=config.A_matrix(),
        node=("v13_seed",0),
    )


def _target(config,provider):
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
    rho0=normalized_grid_density(psi,grid.dx)

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
    rhof=normalized_grid_density(
        exact["psi"][-1],grid.dx
    )
    return grid,psi,rho0,V,rhof


def _dictionary(config):
    return generate_gaussian_dictionary(
        config.q_array(),
        config.p_array(),
        config.A_matrix(),
        config.state,
        cartesian_offsets_2d(
            radius=1.0,
            spacing=0.2,
        ),
        width_scales=(1.0,1.5,2.0,3.0,4.0,6.0),
    )


def _projection_row(build):
    projection=build.projection
    if build.history:
        density_error=build.history[-1].density_error
    else:
        density_error=np.nan

    return {
        "basis_size":len(build.basis),
        "projection_fidelity":float(projection.fidelity),
        "relative_residual":float(projection.relative_residual),
        "initial_density_error":float(density_error),
        "condition_number":float(projection.condition_number),
        "selected_labels":[
            step.selected_label
            for step in build.history
        ],
    }


def _evaluate_reference(
    config,
    provider,
    grid,
    V,
    target_final_density,
    build,
):
    projection=build.projection

    exact_projected=run_exact_2d(
        projection.projected_wavefunction,
        grid.dx,grid.dx,V,
        mass=config.mass,
        dt=0.0025,
        steps=int(round(config.final_time/0.0025)),
        store_every=int(round(config.final_time/0.0025)),
    )
    rho_exact_projected=normalized_grid_density(
        exact_projected["psi"][-1],
        grid.dx,
    )

    gaussian=run_spinor_complete_lvc_gaussians(
        build.basis,
        C0=projection.coefficients,
        provider=provider,
        dt=0.005,
        steps=int(round(config.final_time/0.005)),
        integrator="cayley",
        spawn_action_threshold=1e9,
        max_basis=len(build.basis),
        condition_limit=1e12,
        store_every=max(
            1,
            int(round(config.final_time/0.005))//6,
        ),
    )
    rho=spinor_complete_reduced_density(
        gaussian["final_coefficients"],
        gaussian["final_basis"],
        normalize=True,
    )

    target_pop=density_matrix_populations(
        target_final_density
    )
    pop=density_matrix_populations(rho)

    return {
        "rho":rho,
        "projected_exact_rho":rho_exact_projected,
        "populations":pop,
        "target_populations":target_pop,
        "projected_dynamics_density_error":float(np.linalg.norm(
            rho-rho_exact_projected,ord="fro"
        )),
        "target_density_error":float(np.linalg.norm(
            rho-target_final_density,ord="fro"
        )),
        "target_trace_distance":density_trace_distance(
            rho,target_final_density
        ),
        "target_population_error":float(np.linalg.norm(
            pop-target_pop
        )),
        "purity":density_matrix_purity(rho),
        "target_purity":density_matrix_purity(
            target_final_density
        ),
        "purity_error":float(abs(
            density_matrix_purity(rho)
            -density_matrix_purity(target_final_density)
        )),
        "coherence":complex(rho[0,1]),
        "target_coherence":complex(
            target_final_density[0,1]
        ),
        "coherence_magnitude_error":coherence_magnitude_error(
            rho,target_final_density
        ),
        "coherence_phase_error":coherence_phase_error(
            rho,target_final_density
        ),
        "bloch_vector_error":bloch_vector_error(
            rho,target_final_density
        ),
        "max_norm_drift":float(max(
            abs(r["norm"]-1.0)
            for r in gaussian["records"]
        )),
        "max_condition_number":float(max(
            r["condition_number_nuclear"]
            for r in gaussian["records"]
        )),
    }


def load_v12_context(repository_root):
    path=Path(repository_root)/"results"/"v012_representation_consistent_campaign.json"
    if not path.exists():
        return None

    data=json.loads(path.read_text(encoding="utf-8"))
    r=data["reference_case"]
    return {
        "basis_size":r["n_gaussians"],
        "projection_fidelity":r["initial_projection_fidelity"],
        "relative_residual":r["initial_projection_relative_residual"],
        "initial_density_error":r["initial_density_error"],
        "projected_dynamics_density_error":
            r["projected_dynamics_density_error"],
        "target_density_error":r["target_density_error"],
        "target_population_error":r["target_population_error"],
        "purity_error":r["purity_error"],
        "coherence_phase_error":r["coherence_phase_error"],
        "max_norm_drift":r["max_norm_drift"],
        "max_condition_number":r["max_condition_number"],
        "acceptance":data["acceptance"],
    }


def evaluate_v13_acceptance(
    reference,
    defect_enrichment,
    history,
    thresholds=None,
):
    t=thresholds or V13AcceptanceThresholds()

    predicted=defect_enrichment["predicted_squared_reduction"]
    actual=defect_enrichment["actual_squared_reduction"]
    rel_prediction=abs(actual-predicted)/max(abs(predicted),1e-30)

    monotone=all(
        history[i+1]["relative_residual"]
        < history[i]["relative_residual"]
        for i in range(len(history)-1)
    )

    checks={
        "monotone_residual_refinement":bool(monotone),
        "initial_density_representation":
            reference["initial_density_error"] <= t.max_initial_density_error,
        "projected_dynamics":
            reference["projected_dynamics_density_error"]
            <= t.max_projected_dynamics_density_error,
        "target_full_density":
            reference["target_density_error"]
            <= t.max_target_density_error,
        "target_population":
            reference["target_population_error"]
            <= t.max_target_population_error,
        "coherence_phase":
            reference["coherence_phase_error"] is not None
            and reference["coherence_phase_error"]
            <= t.max_coherence_phase_error,
        "norm":
            reference["max_norm_drift"] <= t.max_norm_drift,
        "conditioning":
            reference["max_condition_number"] <= t.max_condition_number,
        "defect_reduction":
            actual >= t.min_defect_squared_reduction,
        "defect_gain_prediction":
            rel_prediction <= t.max_defect_prediction_relative_error,
    }

    return {
        "passed":bool(all(checks.values())),
        "checks":checks,
        "thresholds":asdict(t),
        "defect_prediction_relative_error":float(rel_prediction),
    }


def run_v013_release_benchmark(
    config=CIPassageConfig(),
    repository_root=None,
):
    provider=AnalyticCI2DFrameProvider(
        nuclear_mass_au=config.mass
    )
    grid,psi,rho0,V,rho_target=_target(
        config,provider
    )
    candidates=_dictionary(config)
    seed=_seed(config)

    prepared=prepare_gaussian_dictionary(
        candidates,
        grid.points,
        grid.dx,
    )

    pure=build_residual_greedy_basis_prepared(
        psi,
        grid.points,
        grid.dx,
        provider,
        [seed],
        prepared,
        max_basis=11,
        top_k_density_screen=1,
        density_screen=False,
        condition_limit=1e5,
    )

    screened=build_residual_greedy_basis_prepared(
        psi,
        grid.points,
        grid.dx,
        provider,
        [seed],
        prepared,
        max_basis=11,
        top_k_density_screen=30,
        density_screen=True,
        condition_limit=1e5,
    )

    reference=_evaluate_reference(
        config,
        provider,
        grid,
        V,
        rho_target,
        screened,
    )
    reference["basis_size"]=len(screened.basis)
    reference["projection_fidelity"]=float(
        screened.projection.fidelity
    )
    reference["relative_residual"]=float(
        screened.projection.relative_residual
    )
    rho_proj0=normalized_grid_density(
        screened.projection.projected_wavefunction,
        grid.dx,
    )
    reference["initial_density_error"]=float(
        np.linalg.norm(rho_proj0-rho0,ord="fro")
    )
    reference["initial_condition_number"]=float(
        screened.projection.condition_number
    )

    defect_result=enrich_basis_from_tdse_defect(
        screened.projection.coefficients,
        screened.basis,
        provider,
        grid,
        candidates,
        condition_limit=1e5,
    )
    if defect_result is None:
        raise RuntimeError(
            "No admissible TDSE-defect enrichment candidate was found."
        )

    defect_enrichment={
        "selected_label":
            defect_result.selected_candidate.label,
        "defect_norm_before":
            defect_result.defect_before.residual_norm,
        "defect_norm_after":
            defect_result.defect_after.residual_norm,
        "relative_defect_before":
            defect_result.defect_before.relative_to_hpsi,
        "relative_defect_after":
            defect_result.defect_after.relative_to_hpsi,
        "predicted_squared_reduction":
            defect_result.score.captured_defect_norm**2,
        "actual_squared_reduction":
            defect_result.actual_squared_defect_reduction,
        "capture_fraction":
            defect_result.score.capture_fraction,
        "expanded_condition_number":
            defect_result.score.expanded_condition_number,
        "zero_coefficient_insertion":True,
    }

    seed_projection=project_grid_wavefunction_to_spinor_complete_basis(
        psi,
        grid.points,
        grid.dx,
        [seed],
        provider,
    )
    seed_rho=normalized_grid_density(
        seed_projection.projected_wavefunction,
        grid.dx,
    )

    ladder=[{
        "basis_size":1,
        "projection_fidelity":float(seed_projection.fidelity),
        "relative_residual":float(seed_projection.relative_residual),
        "density_error":float(np.linalg.norm(
            seed_rho-rho0,ord="fro"
        )),
        "condition_number":float(
            seed_projection.condition_number
        ),
        "selected_label":"seed",
        "predicted_gain":0.0,
        "actual_residual_reduction":0.0,
    }]
    ladder.extend([
        {
            "basis_size":step.basis_size,
            "projection_fidelity":step.projection_fidelity,
            "relative_residual":step.relative_residual,
            "density_error":step.density_error,
            "condition_number":step.condition_number,
            "selected_label":step.selected_label,
            "predicted_gain":step.predicted_gain,
            "actual_residual_reduction":
                step.actual_residual_reduction,
        }
        for step in screened.history
    ])

    pure_summary=_projection_row(pure)
    screened_summary=_projection_row(screened)

    acceptance=evaluate_v13_acceptance(
        reference,
        defect_enrichment,
        ladder,
    )

    context=None
    if repository_root is not None:
        context=load_v12_context(repository_root)

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
        "dictionary":{
            "candidate_count":len(candidates),
            "position_radius":1.0,
            "position_spacing":0.2,
            "width_scales":[1.0,1.5,2.0,3.0,4.0,6.0],
            "momentum_offsets":[[0.0,0.0]],
        },
        "pure_residual_greedy":pure_summary,
        "density_screened_residual":screened_summary,
        "selection_ladder":ladder,
        "reference":reference,
        "defect_enrichment":defect_enrichment,
        "acceptance":acceptance,
        "v12_context":context,
    }


# Imported late to keep the top-level dependency list readable.
from .initial_projection_v12 import (
    project_grid_wavefunction_to_spinor_complete_basis,
)
