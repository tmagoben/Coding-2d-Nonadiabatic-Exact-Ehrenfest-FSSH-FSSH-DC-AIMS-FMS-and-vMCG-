from dataclasses import dataclass, asdict
from pathlib import Path
import json
import time
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
from .sparse_adaptive_dynamics_v16 import (
    SparseAdaptiveSettingsV16,
    run_sparse_cost_aware_lvc_gaussians,
)
from .sparse_pair_matrices_v16 import (
    sparse_reduced_density,
    build_sparse_spinor_lvc_matrices,
    audit_sparse_lvc_matrices_against_dense,
)
from .locality_graph_v16 import (
    LocalityGraphSettings,
    PersistentGaussianLocalityGraph,
)
from .electronic_cost_v16 import (
    GeometryCacheElectronicCostModel,
)
from .local_cost_aware_v16 import (
    estimate_local_sparse_incremental_cost,
)
from .sparse_complexity_v16 import (
    sparse_complexity_model_v16,
)
from .pair_cache_v15 import (
    GaussianPairCache,
    build_cached_spinor_lvc_matrices,
)


@dataclass(frozen=True)
class V16AcceptanceThresholds:
    max_initial_density_error: float = 0.035
    max_projected_dynamics_density_error: float = 0.001
    max_target_density_error: float = 0.035
    max_target_population_error: float = 0.03
    max_coherence_phase_error: float = 0.0035
    max_norm_drift: float = 1e-4
    max_condition_number: float = 5e3

    min_average_graph_sparsity: float = 0.04
    min_pair_factorization_reduction_vs_v15: float = 0.04
    max_final_rho_difference_vs_v15: float = 0.0015

    min_scaling_pair_reduction_n80: float = 0.90
    max_scaling_edge_fraction_n80: float = 0.08
    max_local_edge_scaling_exponent: float = 1.20
    min_dense_pair_scaling_exponent: float = 1.80

    max_final_sparse_S_relative_error: float = 0.01
    max_final_sparse_H_relative_error: float = 0.01

    electronic_cache_demo_requires_lower_cost: bool = True


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
        node=("v16_seed",0),
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


def load_v15_context(repository_root):
    path=Path(repository_root)/"results"/"v015_cost_aware_cache_campaign.json"
    if not path.exists():
        return None

    data=json.loads(
        path.read_text(encoding="utf-8")
    )
    r=data["reference"]
    p=np.asarray(r["populations"],dtype=float)
    coh=complex(*r["coherence"])
    rho=np.array([
        [p[0],coh],
        [np.conj(coh),p[1]],
    ],dtype=complex)

    return {
        "reference":r,
        "final_density_matrix":rho,
        "complexity":
            data["adaptive"]["complexity"],
        "timing_comparison":
            data.get("timing_comparison"),
    }


def _chain_basis(n,spacing=1.5):
    basis=[]
    center=0.5*(n-1)
    for i in range(int(n)):
        basis.append(
            DynamicGraphTBF(
                uid=i,
                state=i%2,
                q=np.array([
                    spacing*(i-center),
                    0.15*np.sin(0.3*i),
                ]),
                p=np.array([
                    0.05*np.cos(0.2*i),
                    0.02*np.sin(0.4*i),
                ]),
                A=np.array([
                    [1.0+0.1*(i%3),0.0],
                    [0.0,0.9+0.05*(i%4)],
                ]),
                node=("chain",i),
            )
        )
    return basis


def run_sparse_scaling_benchmark(
    sizes=(20,40,80),
    *,
    enter_overlap=0.03,
    exit_overlap=0.015,
):
    provider=AnalyticCI2DFrameProvider(
        nuclear_mass_au=20.0
    )
    rows=[]

    for n in sizes:
        basis=_chain_basis(n)
        graph=PersistentGaussianLocalityGraph(
            LocalityGraphSettings(
                enter_overlap=enter_overlap,
                exit_overlap=exit_overlap,
            )
        )

        t0=time.perf_counter()
        update=graph.update(basis)
        graph_seconds=time.perf_counter()-t0

        t0=time.perf_counter()
        mats=build_sparse_spinor_lvc_matrices(
            update,provider
        )
        matrix_seconds=time.perf_counter()-t0

        t0=time.perf_counter()
        dense_cache=GaussianPairCache(basis)
        build_cached_spinor_lvc_matrices(
            dense_cache,provider
        )
        dense_matrix_seconds=time.perf_counter()-t0

        dense_canonical=n*(n+1)//2
        actual_factorizations=(
            update.cache.stats.canonical_solves
        )
        reduction=1.0-(
            actual_factorizations
            /max(dense_canonical,1)
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
            "screened_pairs":
                int(update.screened_pairs),
            "spatial_candidate_pairs":
                int(update.spatial_candidate_pairs),
            "globally_screened_pairs":
                int(update.globally_screened_pairs),
            "exact_pair_checks":
                int(update.exact_pair_checks),
            "pair_factorizations":
                int(actual_factorizations),
            "dense_canonical_pairs":
                int(dense_canonical),
            "pair_reduction_fraction":
                float(reduction),
            "S_nnz":int(mats.S.nnz),
            "H_nnz":int(mats.H.nnz),
            "H_density":float(mats.H_density),
            "graph_seconds":
                float(graph_seconds),
            "matrix_seconds":
                float(matrix_seconds),
            "dense_matrix_seconds":
                float(dense_matrix_seconds),
            "assembly_speedup_vs_dense":
                float(
                    dense_matrix_seconds
                    /max(matrix_seconds,1e-30)
                ),
        })

    return rows


def fit_sparse_scaling_exponents(rows):
    N=np.asarray([r["n_basis"] for r in rows],dtype=float)
    active=np.asarray([r["active_edges"] for r in rows],dtype=float)
    spatial=np.asarray([r["spatial_candidate_pairs"] for r in rows],dtype=float)
    factors=np.asarray([r["pair_factorizations"] for r in rows],dtype=float)
    dense=np.asarray([r["dense_canonical_pairs"] for r in rows],dtype=float)

    def slope(y):
        return float(np.polyfit(
            np.log(N),np.log(y),1
        )[0])

    return {
        "active_edge_exponent":slope(active),
        "spatial_candidate_exponent":slope(spatial),
        "pair_factorization_exponent":slope(factors),
        "dense_canonical_pair_exponent":slope(dense),
        "interpretation":(
            "These fitted exponents apply only to the bounded-locality synthetic "
            "Gaussian chain, not to arbitrary dense Gaussian configurations."
        ),
    }


def electronic_cost_demo():
    basis=[
        DynamicGraphTBF(
            uid=0,state=0,
            q=np.array([0.0,0.0]),
            p=np.zeros(2),
            A=np.eye(2),
            node=("demo",0),
        ),
        DynamicGraphTBF(
            uid=1,state=0,
            q=np.array([0.6,0.0]),
            p=np.zeros(2),
            A=np.eye(2),
            node=("demo",1),
        ),
    ]

    cached=DynamicGraphTBF(
        uid=10,state=1,
        q=np.array([0.04,0.0]),
        p=np.zeros(2),
        A=np.eye(2),
        node=("demo","cached"),
    )
    fresh=DynamicGraphTBF(
        uid=11,state=1,
        q=np.array([1.0,0.0]),
        p=np.zeros(2),
        A=np.eye(2),
        node=("demo","fresh"),
    )

    model=GeometryCacheElectronicCostModel(
        [[0.0,0.0]],
        reuse_radius=0.1,
        cached_cost_units=0.05,
        new_cost_units=2.0,
    )

    common=dict(
        basis=basis,
        active_offdiagonal_edges=1,
        overlap_threshold=0.03,
        horizon_steps=10,
        current_condition=10.0,
        expanded_condition=10.0,
        electronic_cost_model=model,
        electronic_cost_weight=1.0,
    )
    a=estimate_local_sparse_incremental_cost(
        cached,**common
    )
    b=estimate_local_sparse_incremental_cost(
        fresh,**common
    )

    return {
        "cached_geometry":{
            "q":cached.q.tolist(),
            "cost_units":
                a.electronic_cost_units,
            "cache_hit":
                a.electronic_cache_hit,
            "normalized_incremental_cost":
                a.normalized_incremental_cost,
        },
        "new_geometry":{
            "q":fresh.q.tolist(),
            "cost_units":
                b.electronic_cost_units,
            "cache_hit":
                b.electronic_cache_hit,
            "normalized_incremental_cost":
                b.normalized_incremental_cost,
        },
    }


def evaluate_v16_acceptance(
    reference,
    adaptive,
    scaling,
    electronic_demo,
    v15_context=None,
    thresholds=None,
    scaling_fit=None,
    sparse_audit=None,
):
    t=thresholds or V16AcceptanceThresholds()
    c=adaptive["complexity"]

    rho_diff=None
    pair_reduction=None
    if v15_context is not None:
        rho_diff=float(np.linalg.norm(
            reference["final_density_matrix"]
            -v15_context["final_density_matrix"],
            ord="fro",
        ))
        old_pairs=float(
            v15_context["complexity"][
                "propagation_pair_factorizations"
            ]
        )
        pair_reduction=float(
            1.0
            -c["propagation_pair_factorizations"]
            /max(old_pairs,1.0)
        )

    row80=next(
        x for x in scaling
        if x["n_basis"]==80
    )

    demo_ok=(
        electronic_demo["cached_geometry"][
            "normalized_incremental_cost"
        ]
        <
        electronic_demo["new_geometry"][
            "normalized_incremental_cost"
        ]
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
        "graph_is_actually_sparse":
            c["average_sparsity_fraction"]
            >=t.min_average_graph_sparsity,
        "pair_work_reduced_vs_v15":
            pair_reduction is None
            or pair_reduction
            >=t.min_pair_factorization_reduction_vs_v15,
        "sparse_result_close_to_v15":
            rho_diff is None
            or rho_diff
            <=t.max_final_rho_difference_vs_v15,
        "n80_pair_reduction":
            row80["pair_reduction_fraction"]
            >=t.min_scaling_pair_reduction_n80,
        "n80_edge_fraction":
            row80["edge_fraction"]
            <=t.max_scaling_edge_fraction_n80,
        "local_edge_scaling":
            scaling_fit is None
            or scaling_fit["active_edge_exponent"]
            <=t.max_local_edge_scaling_exponent,
        "dense_pair_scaling":
            scaling_fit is None
            or scaling_fit["dense_canonical_pair_exponent"]
            >=t.min_dense_pair_scaling_exponent,
        "final_sparse_S_audit":
            sparse_audit is None
            or sparse_audit["relative_S_frobenius_error"]
            <=t.max_final_sparse_S_relative_error,
        "final_sparse_H_audit":
            sparse_audit is None
            or sparse_audit["relative_H_frobenius_error"]
            <=t.max_final_sparse_H_relative_error,
        "electronic_cache_cost_demo":
            (not t.electronic_cache_demo_requires_lower_cost)
            or demo_ok,
    }

    return {
        "passed":bool(all(checks.values())),
        "checks":checks,
        "thresholds":asdict(t),
        "final_rho_difference_vs_v15":rho_diff,
        "pair_factorization_reduction_vs_v15":
            pair_reduction,
    }


def run_v016_release_benchmark(
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

    control=SparseAdaptiveSettingsV16(
        defect_interval=10,
        enrich_relative_threshold=0.020,
        prune_relative_threshold=0.006,
        minimum_capture_fraction=0.003,
        minimum_local_utility=0.08,
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
        candidate_position_shifts=(0.0,0.06,-0.06),
        candidate_width_scales=(0.75,1.0,1.35),
        locality_enter_overlap=0.03,
        locality_exit_overlap=0.015,
        check_initial_defect=False,
    )

    adaptive=run_sparse_cost_aware_lvc_gaussians(
        build.basis,
        C0=build.projection.coefficients,
        provider=provider,
        grid=defect_grid,
        dt=0.005,
        steps=int(round(config.final_time/0.005)),
        settings=control,
        store_every=10,
    )

    rho=sparse_reduced_density(
        adaptive["final_coefficients"],
        adaptive["final_sparse_matrices"].Snuc,
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
            abs(x["norm"]-1.0)
            for x in adaptive["records"]
        )),
        "max_condition_number":float(max(
            x["condition_number"]
            for x in adaptive["records"]
        )),
        "final_density_matrix":rho,
    }

    v15=None
    if repository_root is not None:
        v15=load_v15_context(repository_root)

    scaling=run_sparse_scaling_benchmark()
    scaling_fit=fit_sparse_scaling_exponents(
        scaling
    )
    demo=electronic_cost_demo()
    sparse_audit=audit_sparse_lvc_matrices_against_dense(
        adaptive["final_basis"],
        provider,
        adaptive["final_sparse_matrices"],
    )

    acceptance=evaluate_v16_acceptance(
        reference,
        {
            "events":adaptive["events"],
            "complexity":adaptive["complexity"],
        },
        scaling,
        demo,
        v15_context=v15,
        scaling_fit=scaling_fit,
        sparse_audit=sparse_audit,
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
            "defect_history":
                adaptive["defect_history"],
            "cost_history":
                adaptive["cost_history"],
            "average_basis_size":
                adaptive["average_basis_size"],
            "endpoint_graph":
                adaptive["endpoint_graph"],
            "midpoint_graph":
                adaptive["midpoint_graph"],
            "complexity":
                adaptive["complexity"],
            "control":asdict(control),
        },
        "sparse_scaling":scaling,
        "sparse_scaling_fit":scaling_fit,
        "final_sparse_matrix_audit":sparse_audit,
        "electronic_cost_demo":demo,
        "v15_context":v15,
        "complexity_model":
            sparse_complexity_model_v16().__dict__.copy(),
        "acceptance":acceptance,
    }
