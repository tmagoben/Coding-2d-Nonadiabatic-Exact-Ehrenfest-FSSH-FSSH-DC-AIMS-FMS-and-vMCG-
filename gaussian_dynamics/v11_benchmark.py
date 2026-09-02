from dataclasses import dataclass, asdict
import gc
import numpy as np

from .benchmark_campaign import (
    CIPassageConfig,
    run_exact_passage,
    run_managed_passage,
)
from .dynamic_graph_aims import DynamicGraphTBF
from .dynamic_gauge_graph import AnalyticCI2DFrameProvider
from .managed_graph_aims_v11 import run_basis_complete_graph_aims
from .electronic_observables import (
    exact_reduced_electronic_density_diabatic,
    reduced_electronic_density_analytic_ci_diabatic,
    density_matrix_populations,
    density_matrix_purity,
    density_matrix_linear_entropy,
    density_matrix_von_neumann_entropy,
)
from .basis_completeness import basis_completeness_report


@dataclass(frozen=True)
class V11AcceptanceThresholds:
    max_population_l2_error: float = 0.05
    max_density_frobenius_error: float = 0.10
    max_purity_error: float = 0.05
    max_norm_drift: float = 1e-2
    max_condition_number: float = 1e6


def _seed(config):
    return DynamicGraphTBF(
        uid=0,
        state=config.state,
        q=config.q_array(),
        p=config.p_array(),
        A=config.A_matrix(),
        node=("seed",0),
    )


def _exact_density(config):
    exact=run_exact_passage(
        config,
        grid_n=64,
        dt=0.0025,
    )
    rho=exact_reduced_electronic_density_diabatic(
        exact["psi_final"],
        exact["dx"],
        exact["dx"],
    )
    rho=rho/np.trace(rho)
    return exact,rho


def _managed_density(run):
    return reduced_electronic_density_analytic_ci_diabatic(
        run["final_coefficients"],
        run["final_basis"],
        normalize=True,
    )


def _metrics(run,rho_exact):
    rho=_managed_density(run)
    p=density_matrix_populations(rho)
    pex=density_matrix_populations(rho_exact)

    return {
        "rho":rho,
        "populations":p,
        "population_l2_error":float(np.linalg.norm(p-pex)),
        "density_frobenius_error":float(np.linalg.norm(rho-rho_exact,ord="fro")),
        "purity":density_matrix_purity(rho),
        "purity_error":float(abs(
            density_matrix_purity(rho)-density_matrix_purity(rho_exact)
        )),
        "linear_entropy":density_matrix_linear_entropy(rho),
        "von_neumann_entropy":density_matrix_von_neumann_entropy(rho),
        "max_norm_drift":float(max(
            abs(float(r["norm"])-1.0) for r in run["records"]
        )),
        "max_condition_number":float(max(
            float(r["condition_number"]) for r in run["records"]
        )),
        "basis_size":len(run["final_basis"]),
    }


def _run_v11(config,max_basis=10,position_shifts=(0.0,0.05,-0.05),width_scales=(0.65,1.0,1.55),children_per_event=2,store_every=120):
    provider=AnalyticCI2DFrameProvider(nuclear_mass_au=config.mass)

    return run_basis_complete_graph_aims(
        [_seed(config)],
        [1.0+0j],
        provider=provider,
        dt=0.005,
        steps=int(round(config.final_time/0.005)),
        spa_order=1,
        spawn_action_threshold=1e-4,
        overlap_block=0.9999,
        child_overlap_block=0.995,
        max_basis=max_basis,
        max_generation=5,
        children_per_event=children_per_event,
        allow_repeated_spawning=True,
        minimum_spawn_separation_steps=4,
        position_shifts=position_shifts,
        width_scales=width_scales,
        momentum_directions=("nac","momentum"),
        condition_limit=1e9,
        eigenvalue_floor=1e-10,
        max_pruning_loss=1e-7,
        store_every=int(store_every),
    )




def run_v011_case(
    config=CIPassageConfig(),
    max_basis=10,
    position_shifts=(0.0,0.05,-0.05),
    width_scales=(0.65,1.0,1.55),
    children_per_event=2,
    store_every=120,
):
    """Public convenience wrapper for one v0.11 strong-CI managed case."""
    return _run_v11(
        config,
        max_basis=max_basis,
        position_shifts=position_shifts,
        width_scales=width_scales,
        children_per_event=children_per_event,
        store_every=store_every,
    )


def _v10_baseline(config):
    return run_managed_passage(
        config=config,
        dt=0.005,
        spa_order=0,
        spawn_action_threshold=2e-4,
        max_basis=4,
        overlap_block=0.9999,
        minimum_spawn_separation_steps=5,
        store_every=20,
    )


def evaluate_v11_acceptance(metrics,thresholds=None):
    t=thresholds or V11AcceptanceThresholds()

    checks={
        "population":metrics["population_l2_error"] <= t.max_population_l2_error,
        "full_density":metrics["density_frobenius_error"] <= t.max_density_frobenius_error,
        "purity":metrics["purity_error"] <= t.max_purity_error,
        "norm":metrics["max_norm_drift"] <= t.max_norm_drift,
        "conditioning":metrics["max_condition_number"] <= t.max_condition_number,
    }

    return {
        "passed":bool(all(checks.values())),
        "checks":checks,
        "thresholds":asdict(t),
    }


def run_v011_release_benchmark(config=CIPassageConfig(), include_ablations=False):
    """Run the deterministic v0.11 release campaign with bounded memory use."""
    exact,rho_exact=_exact_density(config)

    exact_summary={
        "rho":rho_exact,
        "populations":density_matrix_populations(rho_exact),
        "purity":density_matrix_purity(rho_exact),
        "linear_entropy":density_matrix_linear_entropy(rho_exact),
        "von_neumann_entropy":density_matrix_von_neumann_entropy(rho_exact),
    }
    del exact
    gc.collect()

    baseline=_v10_baseline(config)
    baseline_metrics=_metrics(baseline,rho_exact)
    del baseline
    gc.collect()

    basis_rows=[]
    reference_run=None
    reference_metrics=None

    for max_basis in (2,4,6,8,10):
        run=_run_v11(
            config,
            max_basis=max_basis,
            store_every=(20 if max_basis==10 else 120),
        )
        metrics=_metrics(run,rho_exact)

        basis_rows.append({
            "max_basis":int(max_basis),
            **{k:v for k,v in metrics.items() if k!="rho"},
        })

        if max_basis==10:
            reference_run=run
            reference_metrics=metrics
        else:
            del run
        gc.collect()

    # Compact the reference before running expensive ablations.  Dynamic gauge graphs
    # are intentionally not retained in the release JSON.
    completeness=basis_completeness_report(reference_run)
    reference_lineage={
        int(uid):dict(info)
        for uid,info in reference_run["lineage"].items()
    }
    reference_events=list(reference_run["events"])
    reference_settings=dict(reference_run["settings"])
    reference_rho=np.asarray(reference_metrics["rho"]).copy()

    del reference_run
    gc.collect()

    ablations={}

    if include_ablations:
        run=_run_v11(
            config,
            max_basis=10,
            position_shifts=(0.0,),
            store_every=120,
        )
        values=_metrics(run,rho_exact)
        ablations["no_position_optimization"]={
            k:v for k,v in values.items() if k!="rho"
        }
        del run,values
        gc.collect()

        run=_run_v11(
            config,
            max_basis=10,
            width_scales=(1.0,),
            store_every=120,
        )
        values=_metrics(run,rho_exact)
        ablations["fixed_width_only"]={
            k:v for k,v in values.items() if k!="rho"
        }
        del run,values
        gc.collect()

    acceptance=evaluate_v11_acceptance(reference_metrics)

    result={
        "config":{
            "q0":list(config.q0),
            "p0":list(config.p0),
            "A_diag":list(config.A_diag),
            "state":config.state,
            "mass":config.mass,
            "final_time":config.final_time,
            "half_width":config.half_width,
        },
        "exact":exact_summary,
        "v10_baseline":{
            k:v for k,v in baseline_metrics.items() if k!="rho"
        },
        "v11_reference":{
            k:v for k,v in reference_metrics.items() if k!="rho"
        },
        "v11_reference_rho":reference_rho,
        "basis_ladder":basis_rows,
        "ablations":ablations,
        "basis_completeness":completeness,
        "acceptance":acceptance,
        "lineage":reference_lineage,
        "events":reference_events,
        "settings":reference_settings,
    }

    return result

