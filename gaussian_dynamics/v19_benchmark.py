from dataclasses import dataclass, asdict
from pathlib import Path
import importlib.util
import time
import numpy as np

from .analytic_molecular_backend_v19 import (
    AnalyticMolecularLVCBackendV19,
    AnalyticMolecularLVCConfigV19,
    default_diatomic_two_mode_map_v19,
)
from .benchmark_provider_nd import LVC2DGeneralizedProvider
from .molecular_backend import GeneralizedCoordinateProvider
from .molecular_direct_provider_v19 import (
    TrackedMolecularDirectProviderV19,
    BackendEvaluationPolicyV19,
)
from .local_gaussian_nd import LocalAdiabaticTBF
from .molecular_gauge_graph_v19 import (
    build_molecular_centroid_graph_v19,
)
from .molecular_direct_dynamics_v19 import (
    run_molecular_direct_dynamics_v19,
)
from .direct_dynamics_nd import (
    run_backend_spawned_gaussians,
)
from .state_tracking_v19 import (
    scalable_maximum_overlap_assignment_v19,
)


@dataclass(frozen=True)
class V19AcceptanceThresholds:
    max_energy_error: float=1e-11
    max_gradient_error: float=1e-11
    max_nac_error: float=1e-11
    min_raw_scramble_error: float=1e-4

    max_graph_S_difference: float=1e-10
    max_graph_H_difference: float=1e-10
    max_graph_hermiticity_error: float=1e-10

    max_dynamics_coefficient_difference: float=1e-10
    max_dynamics_center_difference: float=1e-10
    max_norm_drift: float=1e-10
    required_spawn_events: int=1

    min_cache_hits: int=10
    max_tracking_ambiguities: int=0
    required_fallback_uses: int=1

    min_large_tracking_states: int=16


def _scan_points():
    return [
        np.array([x,0.35+0.03*np.sin(2.0*x)])
        for x in np.linspace(-0.8,0.8,17)
    ]


def _local_basis():
    A=1.1*np.eye(2)
    return [
        LocalAdiabaticTBF(
            0,np.array([-0.6,0.35]),
            np.array([0.3,0.0]),A
        ),
        LocalAdiabaticTBF(
            1,np.array([0.1,0.40]),
            np.array([0.1,0.1]),A
        ),
        LocalAdiabaticTBF(
            0,np.array([0.7,0.32]),
            np.array([-0.2,0.0]),A
        ),
    ]


def _dynamics_basis():
    return [
        LocalAdiabaticTBF(
            1,
            np.array([-0.45,0.35]),
            np.array([20.0,0.0]),
            0.8*np.eye(2),
        )
    ]


def _provider_scan(provider,reference,points):
    rows=[]
    for q in points:
        p=provider.evaluate(q)
        r=reference.evaluate(q)
        rows.append({
            "q":q.tolist(),
            "energy_error":float(
                np.max(np.abs(p.energies-r.energies))
            ),
            "gradient_error":float(
                np.max(np.abs(
                    p.gradients_q-r.gradients_q
                ))
            ),
            "nac_error":float(
                np.max(np.abs(
                    p.nac_q-r.nac_q
                ))
            ),
            "node_id":
                p.metadata.get("v19_node_id"),
        })
    return rows


def _max_scan(rows,key):
    return float(max(row[key] for row in rows))


def run_v019_release_benchmark():
    gmap=default_diatomic_two_mode_map_v19()

    clean_backend=AnalyticMolecularLVCBackendV19(
        gmap
    )
    scrambled_backend=AnalyticMolecularLVCBackendV19(
        gmap,
        AnalyticMolecularLVCConfigV19(
            scramble_roots=True
        ),
    )

    mass=float(
        clean_backend.generalized_mass_matrix_au[
            0,0
        ]
    )
    reference=LVC2DGeneralizedProvider(
        nuclear_mass_au=mass
    )

    clean=TrackedMolecularDirectProviderV19(
        clean_backend,gmap
    )
    scrambled=TrackedMolecularDirectProviderV19(
        scrambled_backend,gmap
    )

    points=_scan_points()

    # The first point defines the tracked labels, just as the initial molecular
    # geometry would in direct dynamics.
    clean_rows=_provider_scan(
        clean,reference,points
    )
    scrambled_rows=_provider_scan(
        scrambled,reference,points
    )

    # After one explicit reference seed, evaluate the same path in a deliberately
    # nonsequential order to exercise branched/centroid-style tracking.
    shuffled_provider=TrackedMolecularDirectProviderV19(
        AnalyticMolecularLVCBackendV19(
            gmap,
            AnalyticMolecularLVCConfigV19(
                scramble_roots=True
            ),
        ),
        gmap,
    )
    shuffled_provider.evaluate(points[0])
    shuffle_order=[
        8,4,12,2,14,6,10,1,
        15,3,13,5,11,7,9,16,
    ]
    shuffled_rows=_provider_scan(
        shuffled_provider,
        reference,
        [points[i] for i in shuffle_order],
    )

    # Repeat in reverse order to exercise exact-geometry cache reuse.
    for q in reversed(points):
        scrambled.evaluate(q)

    # Raw untracked scrambled provider: demonstrates why overlap tracking is needed.
    raw_scrambled=GeneralizedCoordinateProvider(
        AnalyticMolecularLVCBackendV19(
            gmap,
            AnalyticMolecularLVCConfigV19(
                scramble_roots=True
            ),
        ),
        gmap,
    )
    raw_errors=[]
    for q in points:
        p=raw_scrambled.evaluate(q)
        r=reference.evaluate(q)
        raw_errors.append(float(
            np.max(np.abs(
                p.energies-r.energies
            ))
        ))

    # Center-centroid gauge graph.
    clean_graph=build_molecular_centroid_graph_v19(
        _local_basis(),clean
    )
    scrambled_graph=build_molecular_centroid_graph_v19(
        _local_basis(),scrambled
    )
    Sc,Hc=clean_graph.matrices()
    Ss,Hs=scrambled_graph.matrices()

    graph_result={
        "basis_size":3,
        "nodes":len(clean_graph.registry.graph.nodes),
        "edges":len(clean_graph.registry.graph.edges()),
        "S_difference":float(
            np.linalg.norm(Sc-Ss,ord="fro")
        ),
        "H_difference":float(
            np.linalg.norm(Hc-Hs,ord="fro")
        ),
        "S_hermiticity_error":float(
            np.linalg.norm(
                Sc-Sc.conj().T,ord="fro"
            )
        ),
        "H_hermiticity_error":float(
            np.linalg.norm(
                Hc-Hc.conj().T,ord="fro"
            )
        ),
        "condition_number":float(
            np.linalg.cond(Sc)
        ),
    }

    # Direct dynamics: analytic generalized provider, clean molecular provider, and
    # deliberately scrambled+tracked molecular provider must agree.
    dyn_kwargs=dict(
        initial_basis=_dynamics_basis(),
        C0=np.array([1.0+0.0j]),
        dt=0.02,
        steps=50,
        spawn_threshold=0.001,
        overlap_block=0.999999,
        max_basis=2,
        store_every=10,
    )
    analytic_dyn=run_backend_spawned_gaussians(
        provider=reference,
        **dyn_kwargs,
    )
    molecular_dyn=run_molecular_direct_dynamics_v19(
        provider=TrackedMolecularDirectProviderV19(
            AnalyticMolecularLVCBackendV19(gmap),
            gmap,
        ),
        **dyn_kwargs,
    )
    scrambled_dyn_provider=TrackedMolecularDirectProviderV19(
        AnalyticMolecularLVCBackendV19(
            gmap,
            AnalyticMolecularLVCConfigV19(
                scramble_roots=True
            ),
        ),
        gmap,
    )
    scrambled_dyn=run_molecular_direct_dynamics_v19(
        provider=scrambled_dyn_provider,
        **dyn_kwargs,
    )

    coeff_diff=float(np.linalg.norm(
        analytic_dyn["final_coefficients"]
        -scrambled_dyn["final_coefficients"]
    ))
    center_diff=float(max(
        np.linalg.norm(a.q-b.q)
        for a,b in zip(
            analytic_dyn["final_basis"],
            scrambled_dyn["final_basis"],
        )
    ))
    momentum_diff=float(max(
        np.linalg.norm(a.p-b.p)
        for a,b in zip(
            analytic_dyn["final_basis"],
            scrambled_dyn["final_basis"],
        )
    ))
    norm_drift=float(max(
        np.max(np.abs(
            analytic_dyn["norm"]-1.0
        )),
        np.max(np.abs(
            scrambled_dyn["norm"]-1.0
        )),
    ))

    dynamics_result={
        "analytic_spawn_events":
            analytic_dyn["events"],
        "molecular_spawn_events":
            molecular_dyn["events"],
        "scrambled_spawn_events":
            scrambled_dyn["events"],
        "coefficient_difference":
            coeff_diff,
        "center_difference":
            center_diff,
        "momentum_difference":
            momentum_diff,
        "maximum_norm_drift":
            norm_drift,
        "final_coefficients_analytic":
            analytic_dyn["final_coefficients"],
        "final_coefficients_scrambled":
            scrambled_dyn["final_coefficients"],
        "scrambled_provider_diagnostics":
            scrambled_dyn_provider.diagnostics_dict(),
        "final_graph_audit":
            scrambled_dyn[
                "molecular_centroid_graph_audit"
            ],
    }

    # Failure/fallback contract.
    fallback_provider=TrackedMolecularDirectProviderV19(
        AnalyticMolecularLVCBackendV19(
            gmap,
            AnalyticMolecularLVCConfigV19(
                fail_if_q0_greater_than=0.10
            ),
        ),
        gmap,
        failure=BackendEvaluationPolicyV19(
            failure_policy="nearest_cache",
            max_fallback_distance=0.05,
        ),
    )
    fallback_provider.evaluate(
        np.array([0.08,0.30])
    )
    fallback_point=fallback_provider.evaluate(
        np.array([0.12,0.30])
    )

    fallback_result={
        "metadata":
            fallback_point.metadata,
        "diagnostics":
            fallback_provider.diagnostics_dict(),
    }

    # Provider cost interface.
    cost_provider=TrackedMolecularDirectProviderV19(
        AnalyticMolecularLVCBackendV19(gmap),
        gmap,
    )
    q0=np.array([-0.40,0.30])
    cost_provider.evaluate(q0)
    cost_result={
        "cached":
            cost_provider.cost_estimate(q0),
        "nearby":
            cost_provider.cost_estimate(
                q0+np.array([0.01,0.0])
            ),
        "new":
            cost_provider.cost_estimate(
                q0+np.array([0.30,0.0])
            ),
    }

    # Scalable assignment diagnostic.
    rng=np.random.default_rng(190019)
    nlarge=16
    X=rng.normal(size=(nlarge,nlarge))
    Q,_=np.linalg.qr(X)
    t0=time.perf_counter()
    large_assignment=(
        scalable_maximum_overlap_assignment_v19(
            Q,
            minimum_overlap=0.0,
            minimum_score_margin=0.0,
            real_gauge=True,
        )
    )
    assignment_seconds=time.perf_counter()-t0

    tracking_scaling={
        "nstate":nlarge,
        "assignment_seconds":
            float(assignment_seconds),
        "best_score":
            float(large_assignment.best_score),
        "second_best_score":
            float(
                large_assignment.second_best_score
            ),
        "permutation_is_valid":
            sorted(
                large_assignment.permutation.tolist()
            )==list(range(nlarge)),
        "complexity":
            "O(n^4) including exact second-best margin; best matching itself O(n^3)",
        "legacy_complexity":
            "O(n!) exhaustive permutations",
    }

    pyscf_available=(
        importlib.util.find_spec("pyscf")
        is not None
    )

    max_energy=max(
        _max_scan(clean_rows,"energy_error"),
        _max_scan(scrambled_rows,"energy_error"),
    )
    max_grad=max(
        _max_scan(clean_rows,"gradient_error"),
        _max_scan(scrambled_rows,"gradient_error"),
    )
    max_nac=max(
        _max_scan(clean_rows,"nac_error"),
        _max_scan(scrambled_rows,"nac_error"),
    )
    max_shuffled_energy=_max_scan(
        shuffled_rows,"energy_error"
    )
    max_shuffled_nac=_max_scan(
        shuffled_rows,"nac_error"
    )

    thresholds=V19AcceptanceThresholds()
    checks={
        "molecular_energy_projection":
            max_energy<=thresholds.max_energy_error,
        "molecular_gradient_projection":
            max_grad<=thresholds.max_gradient_error,
        "molecular_nac_projection":
            max_nac<=thresholds.max_nac_error,
        "raw_scramble_is_nontrivial":
            max(raw_errors)
            >=thresholds.min_raw_scramble_error,

        "graph_S_gauge_invariance":
            graph_result["S_difference"]
            <=thresholds.max_graph_S_difference,
        "graph_H_gauge_invariance":
            graph_result["H_difference"]
            <=thresholds.max_graph_H_difference,
        "graph_hermiticity":
            max(
                graph_result["S_hermiticity_error"],
                graph_result["H_hermiticity_error"],
            )<=thresholds.max_graph_hermiticity_error,

        "dynamics_coefficients":
            coeff_diff
            <=thresholds.max_dynamics_coefficient_difference,
        "dynamics_centers":
            center_diff
            <=thresholds.max_dynamics_center_difference,
        "dynamics_norm":
            norm_drift
            <=thresholds.max_norm_drift,
        "spawn_event":
            len(scrambled_dyn["events"])
            ==thresholds.required_spawn_events,

        "cache_reuse":
            scrambled.diagnostics.cache_hits
            >=thresholds.min_cache_hits,
        "tracking_unambiguous":
            scrambled.diagnostics.tracking_ambiguities
            <=thresholds.max_tracking_ambiguities,
        "order_tolerant_tracking":
            max_shuffled_energy
            <=thresholds.max_energy_error
            and max_shuffled_nac
            <=thresholds.max_nac_error
            and shuffled_provider.diagnostics.tracking_ambiguities==0,
        "failure_fallback":
            fallback_provider.diagnostics.fallback_uses
            ==thresholds.required_fallback_uses,

        "large_state_assignment":
            nlarge
            >=thresholds.min_large_tracking_states
            and tracking_scaling[
                "permutation_is_valid"
            ],
    }

    return {
        "model":{
            "geometry_map":
                "synthetic H2 two-mode Cartesian embedding",
            "masses_amu":
                list(clean_backend.config.masses_amu),
            "mass_matrix_q_au":
                clean_backend.generalized_mass_matrix_au,
            "nstate":2,
            "nq":2,
            "note":(
                "This is a deterministic molecular-style validation backend, not "
                "an ab-initio molecular calculation."
            ),
        },
        "provider_scan":{
            "points":len(points),
            "clean_rows":clean_rows,
            "scrambled_tracked_rows":
                scrambled_rows,
            "raw_scrambled_max_energy_error":
                float(max(raw_errors)),
            "maximum_tracked_energy_error":
                max_energy,
            "maximum_tracked_gradient_error":
                max_grad,
            "maximum_tracked_nac_error":
                max_nac,
            "shuffled_after_seed_max_energy_error":
                float(max_shuffled_energy),
            "shuffled_after_seed_max_nac_error":
                float(max_shuffled_nac),
            "shuffled_after_seed_diagnostics":
                shuffled_provider.diagnostics_dict(),
            "scrambled_provider_diagnostics":
                scrambled.diagnostics_dict(),
        },
        "gauge_graph":graph_result,
        "direct_dynamics":dynamics_result,
        "failure_fallback":fallback_result,
        "provider_cost":cost_result,
        "state_tracking_scaling":
            tracking_scaling,
        "pyscf":{
            "installed_in_build_environment":
                bool(pyscf_available),
            "runtime_validated":
                False,
            "bridge":
                "PySCFRawSnapshotBackendV19 + casscf_state_overlap_matrix",
        },
        "acceptance":{
            "passed":bool(all(checks.values())),
            "checks":checks,
            "thresholds":asdict(thresholds),
        },
    }
