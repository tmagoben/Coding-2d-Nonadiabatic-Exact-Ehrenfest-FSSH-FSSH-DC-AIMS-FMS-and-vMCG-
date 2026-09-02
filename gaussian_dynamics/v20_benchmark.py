from dataclasses import dataclass, asdict
import importlib.util
import math
import time
import numpy as np

from .analytic_molecular_backend_v19 import (
    AnalyticMolecularLVCBackendV19,
    AnalyticMolecularLVCConfigV19,
    default_diatomic_two_mode_map_v19,
)
from .indexed_molecular_provider_v20 import (
    IndexedTrackedMolecularDirectProviderV20,
)
from .sparse_molecular_matrices_v20 import (
    SparseMolecularTBFV20,
    MolecularSparseSettingsV20,
    SparseMolecularEdgeGraphV20,
    build_sparse_molecular_matrices_v20,
    build_dense_molecular_reference_v20,
    relative_frobenius_error,
)
from .sparse_molecular_dynamics_v20 import (
    SparseMolecularDynamicsSettingsV20,
    run_sparse_molecular_dynamics_v20,
    run_dense_molecular_reference_dynamics_v20,
)


@dataclass(frozen=True)
class V20AcceptanceThresholds:
    max_metric_coefficient_error: float=1e-3
    max_center_error: float=1e-12
    max_norm_drift: float=1e-10

    max_final_S_error: float=0.005
    max_final_H_error: float=0.003
    max_final_T_error: float=0.020
    min_canonical_sparsity: float=0.75
    min_backend_miss_reduction: float=0.40
    max_sampled_audit_failures: int=0

    max_finest_threshold_S_error: float=0.005
    max_finest_threshold_H_error: float=0.003
    max_finest_threshold_T_error: float=0.020

    max_zero_budget_S_error: float=1e-5
    max_zero_budget_H_error: float=1e-5
    max_zero_budget_T_error: float=1e-4

    max_active_edge_exponent: float=1.20
    max_exact_pair_check_exponent: float=1.10
    max_provider_miss_exponent: float=1.10
    min_dense_pair_exponent: float=1.90
    min_n160_pair_check_reduction: float=0.95

    require_controller_relaxation: bool=True


def _provider(
    *,
    scramble=True,
    rebuild_batch=16,
):
    gmap=default_diatomic_two_mode_map_v19()
    return IndexedTrackedMolecularDirectProviderV20(
        AnalyticMolecularLVCBackendV19(
            gmap,
            AnalyticMolecularLVCConfigV19(
                scramble_roots=bool(scramble)
            ),
        ),
        gmap,
        rebuild_batch=int(rebuild_batch),
    )


def _canonical_basis(n=20):
    n=int(n)
    spacing=0.8
    A=2.5*np.eye(2)
    x0=-0.5*spacing*(n-1)
    out=[]
    for i in range(n):
        p0=0.3*np.exp(
            -((i-(n-1)/2)/3.0)**2
        )
        out.append(
            SparseMolecularTBFV20(
                uid=i,
                state=i%2,
                q=np.array([
                    x0+i*spacing,0.35
                ]),
                p=np.array([
                    p0,0.01*(-1)**i
                ]),
                A=A,
            )
        )
    return out


def _canonical_C0(n=20):
    idx=np.arange(int(n),dtype=float)
    center=0.5*(int(n)-1)
    return (
        np.exp(
            -0.5*((idx-center)/2.0)**2
        )
        *np.exp(0.07j*idx)
    ).astype(complex)


def _canonical_graph_settings():
    return MolecularSparseSettingsV20(
        enter_score=0.030,
        exit_score=0.015,
        search_overlap_floor=1e-5,
        local_omitted_score_l2_budget=0.010,
    )


def _canonical_dynamics_settings():
    return SparseMolecularDynamicsSettingsV20(
        graph=_canonical_graph_settings(),
        sampled_audit_interval=5,
        sampled_audit_priority_pairs=4,
        sampled_audit_random_pairs=4,
        dense_sentinel_S_limit=0.01,
        dense_sentinel_H_limit=0.01,
        dense_sentinel_T_limit=0.03,
    )


def _metric_coefficient_error(
    sparse_coeff,
    dense_coeff,
    dense_S,
):
    cs=np.asarray(sparse_coeff,complex)
    cd=np.asarray(dense_coeff,complex)
    S=dense_S
    inner=np.vdot(cd,S@cs)
    phase=(
        1.0+0.0j
        if abs(inner)<1e-30
        else np.exp(-1j*np.angle(inner))
    )
    diff=phase*cs-cd
    num=float(np.real(
        np.vdot(diff,S@diff)
    ))
    den=float(np.real(
        np.vdot(cd,S@cd)
    ))
    return float(
        math.sqrt(max(num,0.0))
        /math.sqrt(max(den,1e-30))
    )


def _snapshot_basis(n=20):
    n=int(n)
    spacing=0.8
    A=2.5*np.eye(2)
    x0=-0.5*spacing*(n-1)
    return [
        SparseMolecularTBFV20(
            i,i%2,
            np.array([x0+i*spacing,0.35]),
            np.array([0.1*(-1)**i,0.0]),
            A,
        )
        for i in range(n)
    ]


def _dense_snapshot_reference(basis,dt):
    settings=MolecularSparseSettingsV20(
        enter_score=1e-14,
        exit_score=1e-14,
        search_overlap_floor=1e-14,
        local_omitted_score_l2_budget=0.0,
        use_kdtree=False,
    )
    return build_dense_molecular_reference_v20(
        basis,_provider(),dt,settings
    )


def _matrix_row(
    basis,
    dense,
    *,
    enter_score,
    budget,
    dt=0.005,
):
    settings=MolecularSparseSettingsV20(
        enter_score=float(enter_score),
        exit_score=0.5*float(enter_score),
        search_overlap_floor=1e-6,
        local_omitted_score_l2_budget=
            float(budget),
    )
    provider=_provider()
    update=SparseMolecularEdgeGraphV20(
        provider,dt,settings
    ).update(basis)
    mats=build_sparse_molecular_matrices_v20(
        basis,update
    )
    return {
        "enter_score":float(enter_score),
        "budget":float(budget),
        "active_edges":
            int(len(update.active_edges)),
        "exact_pair_checks":
            int(update.exact_pair_checks),
        "budget_promoted_edges":
            int(update.budget_promoted_edges),
        "omitted_score_l2":
            float(update.omitted_score_l2),
        "S_error":relative_frobenius_error(
            mats.S.toarray(),dense["S"]
        ),
        "H_error":relative_frobenius_error(
            mats.H.toarray(),dense["H"]
        ),
        "T_error":relative_frobenius_error(
            mats.T_seed.toarray(),
            dense["T_seed"],
        ),
    }


def _nonincreasing(values,tol=1e-14):
    return all(
        b<=a+tol
        for a,b in zip(values[:-1],values[1:])
    )


def _scaling_basis(n):
    n=int(n)
    spacing=2.0
    A=1.4*np.eye(2)
    center=0.5*(n-1)
    out=[]
    for i in range(n):
        x=(
            (i-center)*spacing
            +0.07*np.sin(1.7*i)
        )
        y=0.35+0.025*np.cos(0.91*i)
        out.append(
            SparseMolecularTBFV20(
                i,i%2,
                np.array([x,y]),
                np.array([
                    0.05*(-1)**i,
                    0.01*np.sin(i),
                ]),
                A,
            )
        )
    return out


def _fit_exponent(rows,key):
    x=np.log(np.asarray([
        row["n_basis"] for row in rows
    ],float))
    y=np.log(np.asarray([
        max(float(row[key]),1.0)
        for row in rows
    ],float))
    return float(np.polyfit(x,y,1)[0])


def _scaling_campaign():
    rows=[]
    settings=MolecularSparseSettingsV20(
        enter_score=0.030,
        exit_score=0.015,
        search_overlap_floor=1e-4,
        local_omitted_score_l2_budget=0.010,
    )
    for n in (20,40,80,160):
        provider=_provider(
            rebuild_batch=32
        )
        basis=_scaling_basis(n)
        t0=time.perf_counter()
        update=SparseMolecularEdgeGraphV20(
            provider,0.005,settings
        ).update(basis)
        elapsed=time.perf_counter()-t0
        diag=provider.diagnostics_dict()
        total=int(update.total_offdiagonal_pairs)
        rows.append({
            "n_basis":int(n),
            "active_edges":
                int(len(update.active_edges)),
            "exact_pair_checks":
                int(update.exact_pair_checks),
            "dense_offdiagonal_pairs":total,
            "pair_check_reduction_fraction":
                float(
                    1.0-update.exact_pair_checks
                    /max(total,1)
                ),
            "matrix_sparsity_fraction":
                float(update.sparsity_fraction),
            "provider_cache_misses":
                int(diag["cache_misses"]),
            "provider_cache_hits":
                int(diag["cache_hits"]),
            "index_rebuilds":
                int(
                    diag["spatial_index"][
                        "rebuilds"
                    ]
                ),
            "index_buffer_distance_checks":
                int(
                    diag["spatial_index"][
                        "buffer_distance_checks"
                    ]
                ),
            "wall_seconds":
                float(elapsed),
        })
    fit={
        "active_edge_exponent":
            _fit_exponent(
                rows,"active_edges"
            ),
        "exact_pair_check_exponent":
            _fit_exponent(
                rows,"exact_pair_checks"
            ),
        "provider_miss_exponent":
            _fit_exponent(
                rows,"provider_cache_misses"
            ),
        "dense_pair_exponent":
            _fit_exponent(
                rows,"dense_offdiagonal_pairs"
            ),
    }
    return rows,fit


def _controller_demo():
    basis=_canonical_basis(2)
    aggressive=MolecularSparseSettingsV20(
        enter_score=0.030,
        exit_score=0.015,
        search_overlap_floor=0.90,
        local_omitted_score_l2_budget=0.010,
    )
    settings=SparseMolecularDynamicsSettingsV20(
        graph=aggressive,
        sampled_audit_interval=1,
        sampled_audit_priority_pairs=2,
        sampled_audit_random_pairs=0,
        sampled_audit_search_factor=0.1,
        sampled_audit_relaxation_factor=0.5,
        max_sampled_audit_relaxations=3,
        dense_sentinel_S_limit=1.0,
        dense_sentinel_H_limit=1.0,
        dense_sentinel_T_limit=1.0,
    )
    out=run_sparse_molecular_dynamics_v20(
        basis,
        np.array([1.0+0j,0.05+0j]),
        _provider(),
        dt=0.001,
        steps=1,
        settings=settings,
        store_every=1,
        sentinel_provider_factory=_provider,
    )
    relax=[
        e for e in out["events"]
        if e["kind"]
        =="sampled_molecular_search_relaxation"
    ]
    return {
        "relaxation_events":relax,
        "audit_history":
            out["sampled_audits"],
        "final_audit_passed":bool(
            out["sampled_audits"][-1][
                "passed"
            ]
        ),
    }


def run_v020_release_benchmark():
    # Canonical sparse-vs-dense moving-basis propagation.
    basis=_canonical_basis(20)
    C0=_canonical_C0(20)

    ps=_provider()
    t0=time.perf_counter()
    sparse_out=run_sparse_molecular_dynamics_v20(
        basis,C0,ps,
        dt=0.002,
        steps=20,
        settings=_canonical_dynamics_settings(),
        store_every=5,
        sentinel_provider_factory=_provider,
    )
    sparse_wall=time.perf_counter()-t0

    pd=_provider()
    t0=time.perf_counter()
    dense_out=(
        run_dense_molecular_reference_dynamics_v20(
            basis,C0,pd,
            dt=0.002,
            steps=20,
            store_every=5,
        )
    )
    dense_wall=time.perf_counter()-t0

    metric_error=_metric_coefficient_error(
        sparse_out["final_coefficients"],
        dense_out["final_coefficients"],
        dense_out["final_S"],
    )
    center_error=float(max(
        np.linalg.norm(a.q-b.q)
        for a,b in zip(
            sparse_out["final_basis"],
            dense_out["final_basis"],
        )
    ))
    momentum_error=float(max(
        np.linalg.norm(a.p-b.p)
        for a,b in zip(
            sparse_out["final_basis"],
            dense_out["final_basis"],
        )
    ))
    norm_drift=float(max(
        abs(row["norm"]-1.0)
        for row in sparse_out["records"]
    ))
    avg_sparsity=float(np.mean([
        row["sparsity_fraction"]
        for row in sparse_out["records"]
    ]))
    sampled_failures=int(sum(
        not row["passed"]
        for row in sparse_out[
            "sampled_audits"
        ]
    ))
    sparse_diag=sparse_out[
        "provider_diagnostics"
    ]
    dense_diag=dense_out[
        "provider_diagnostics"
    ]
    backend_miss_reduction=float(
        1.0
        -sparse_diag["cache_misses"]
        /max(dense_diag["cache_misses"],1)
    )

    canonical={
        "metric_coefficient_error":
            metric_error,
        "center_error":center_error,
        "momentum_error":momentum_error,
        "maximum_norm_drift":norm_drift,
        "average_sparsity_fraction":
            avg_sparsity,
        "final_active_edges":
            int(
                sparse_out["records"][-1][
                    "active_edges"
                ]
            ),
        "total_offdiagonal_pairs":
            int(
                sparse_out["records"][-1][
                    "total_pairs"
                ]
            ),
        "graph_total_exact_pair_checks":
            int(
                sparse_out[
                    "graph_total_exact_pair_checks"
                ]
            ),
        "sampled_audit_failures":
            sampled_failures,
        "sampled_audits":
            sparse_out["sampled_audits"],
        "sentinels":
            sparse_out["sentinels"],
        "sparse_provider":{
            "evaluate_calls":
                int(
                    sparse_diag[
                        "evaluate_calls"
                    ]
                ),
            "cache_hits":
                int(sparse_diag["cache_hits"]),
            "cache_misses":
                int(
                    sparse_diag[
                        "cache_misses"
                    ]
                ),
            "spatial_index":
                sparse_diag[
                    "spatial_index"
                ],
        },
        "dense_provider":{
            "evaluate_calls":
                int(
                    dense_diag[
                        "evaluate_calls"
                    ]
                ),
            "cache_hits":
                int(dense_diag["cache_hits"]),
            "cache_misses":
                int(
                    dense_diag[
                        "cache_misses"
                    ]
                ),
        },
        "backend_miss_reduction_fraction":
            backend_miss_reduction,
        "sparse_wall_seconds":
            float(sparse_wall),
        "dense_wall_seconds":
            float(dense_wall),
        "wall_speedup":
            float(
                dense_wall/max(sparse_wall,1e-30)
            ),
        "records":
            sparse_out["records"],
    }

    # Molecular sparse matrix convergence coordinates.
    snapshot_basis=_snapshot_basis(20)
    dense_snapshot=_dense_snapshot_reference(
        snapshot_basis,0.005
    )

    threshold_rows=[
        _matrix_row(
            snapshot_basis,dense_snapshot,
            enter_score=threshold,
            budget=1e9,
        )
        for threshold in (
            0.12,0.08,0.05,
            0.03,0.02,0.01,0.005,
        )
    ]
    budget_rows=[
        _matrix_row(
            snapshot_basis,dense_snapshot,
            enter_score=0.03,
            budget=budget,
        )
        for budget in (
            1e9,0.05,0.02,
            0.01,0.005,0.0,
        )
    ]
    convergence={
        "threshold_rows":
            threshold_rows,
        "budget_rows":
            budget_rows,
        "threshold_monotone":{
            "S":_nonincreasing([
                r["S_error"]
                for r in threshold_rows
            ]),
            "H":_nonincreasing([
                r["H_error"]
                for r in threshold_rows
            ]),
            "T":_nonincreasing([
                r["T_error"]
                for r in threshold_rows
            ]),
        },
        "budget_monotone":{
            "S":_nonincreasing([
                r["S_error"]
                for r in budget_rows
            ]),
            "H":_nonincreasing([
                r["H_error"]
                for r in budget_rows
            ]),
            "T":_nonincreasing([
                r["T_error"]
                for r in budget_rows
            ]),
        },
    }

    scaling,scaling_fit=_scaling_campaign()
    controller=_controller_demo()

    thresholds=V20AcceptanceThresholds()
    final_sentinel=canonical[
        "sentinels"
    ]["final"]
    finest=threshold_rows[-1]
    zero_budget=budget_rows[-1]
    n160=scaling[-1]

    checks={
        "canonical_coefficients":
            metric_error
            <=thresholds.max_metric_coefficient_error,
        "canonical_centers":
            center_error
            <=thresholds.max_center_error,
        "canonical_norm":
            norm_drift
            <=thresholds.max_norm_drift,

        "final_S_sentinel":
            final_sentinel[
                "relative_S_frobenius_error"
            ]<=thresholds.max_final_S_error,
        "final_H_sentinel":
            final_sentinel[
                "relative_H_frobenius_error"
            ]<=thresholds.max_final_H_error,
        "final_T_sentinel":
            final_sentinel[
                "relative_Tseed_frobenius_error"
            ]<=thresholds.max_final_T_error,
        "canonical_sparsity":
            avg_sparsity
            >=thresholds.min_canonical_sparsity,
        "backend_miss_reduction":
            backend_miss_reduction
            >=thresholds.min_backend_miss_reduction,
        "sampled_audits":
            sampled_failures
            <=thresholds.max_sampled_audit_failures,

        "threshold_monotone":
            all(
                convergence[
                    "threshold_monotone"
                ].values()
            ),
        "threshold_finest":
            finest["S_error"]
            <=thresholds.max_finest_threshold_S_error
            and finest["H_error"]
            <=thresholds.max_finest_threshold_H_error
            and finest["T_error"]
            <=thresholds.max_finest_threshold_T_error,

        "budget_monotone":
            all(
                convergence[
                    "budget_monotone"
                ].values()
            ),
        "zero_budget_limit":
            zero_budget["S_error"]
            <=thresholds.max_zero_budget_S_error
            and zero_budget["H_error"]
            <=thresholds.max_zero_budget_H_error
            and zero_budget["T_error"]
            <=thresholds.max_zero_budget_T_error,

        "active_edge_scaling":
            scaling_fit[
                "active_edge_exponent"
            ]<=thresholds.max_active_edge_exponent,
        "pair_check_scaling":
            scaling_fit[
                "exact_pair_check_exponent"
            ]<=thresholds.max_exact_pair_check_exponent,
        "provider_miss_scaling":
            scaling_fit[
                "provider_miss_exponent"
            ]<=thresholds.max_provider_miss_exponent,
        "dense_pair_scaling":
            scaling_fit[
                "dense_pair_exponent"
            ]>=thresholds.min_dense_pair_exponent,
        "n160_pair_reduction":
            n160[
                "pair_check_reduction_fraction"
            ]>=thresholds.min_n160_pair_check_reduction,

        "online_search_controller":
            (
                (not thresholds.require_controller_relaxation)
                or len(
                    controller[
                        "relaxation_events"
                    ]
                )>=1
            )
            and controller[
                "final_audit_passed"
            ],
    }

    return {
        "canonical":canonical,
        "molecular_sparse_convergence":
            convergence,
        "scaling":scaling,
        "scaling_fit":scaling_fit,
        "controller_demo":controller,
        "pyscf":{
            "installed_in_build_environment":
                bool(
                    importlib.util.find_spec(
                        "pyscf"
                    )
                    is not None
                ),
            "runtime_validated":False,
            "note":(
                "v0.20 sparse molecular machinery is validated with the "
                "deterministic Cartesian LVC backend. Real PySCF runtime "
                "validation remains environment dependent."
            ),
        },
        "acceptance":{
            "passed":bool(all(checks.values())),
            "checks":checks,
            "thresholds":asdict(thresholds),
        },
    }
