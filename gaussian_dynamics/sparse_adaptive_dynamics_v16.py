from dataclasses import dataclass, asdict
import numpy as np
from scipy.sparse.linalg import spsolve

from .dynamic_gauge_graph import AnalyticCI2DFrameProvider
from .dynamic_graph_aims import DynamicGraphTBF, _verlet, _kinematics
from .spinor_complete_dynamics_v12 import (
    initialize_spinor_complete_coefficients,
)
from .spinor_complete_lvc_v12 import (
    coefficients_matrix,
    flatten_coefficients,
)
from .tdse_defect_v13 import (
    TDSEDefect,
    reconstruct_spinor_complete_wavefunction,
    reconstruct_spinor_complete_time_derivative,
    apply_lvc_grid_hamiltonian,
)
from .gaussian_nd import gaussian_nd
from .residual_pruning_v14 import (
    prune_low_loss_gaussian_pair,
)
from .defect_candidates_v15 import (
    generate_energy_conserving_defect_candidates_v15,
    rank_dynamic_defect_candidates_cached,
)
from .locality_graph_v16 import (
    LocalityGraphSettings,
    PersistentGaussianLocalityGraph,
)
from .sparse_pair_matrices_v16 import (
    build_sparse_spinor_lvc_matrices,
    build_sparse_spinor_time_matrix,
    sparse_metric_compatible_connection,
    sparse_moving_basis_midpoint_cayley_step,
    sparse_generalized_norm,
    sparse_reduced_density,
)
from .local_cost_aware_v16 import (
    rank_local_sparse_candidates,
)
from .electronic_cost_v16 import (
    UniformElectronicCostModel,
    ElectronicCostEstimate,
)
from .sparse_complexity_v16 import (
    SparseComplexityLedgerV16,
)


@dataclass(frozen=True)
class SparseAdaptiveSettingsV16:
    defect_interval: int = 10
    enrich_relative_threshold: float = 0.020
    prune_relative_threshold: float = 0.006

    minimum_capture_fraction: float = 0.003
    minimum_local_utility: float = 0.08
    condition_penalty_weight: float = 0.15
    electronic_cost_weight: float = 1.0
    cost_horizon_steps: int = 10
    residual_shortlist: int = 8

    min_basis: int = 8
    max_basis: int = 11
    minimum_adaptation_separation_steps: int = 10
    minimum_prune_age_steps: int = 20
    prune_patience_checks: int = 2

    max_prune_fractional_loss: float = 5e-7
    max_replacement_prune_fractional_loss: float = 5e-7
    emergency_prune_fractional_loss: float = 1e-4

    condition_limit: float = 1e5
    hard_condition_limit: float = 5e6
    orthogonal_norm_floor: float = 1e-8

    candidate_position_shifts: tuple = (0.0,0.06,-0.06)
    candidate_width_scales: tuple = (0.75,1.0,1.35)
    candidate_momentum_directions: tuple = ("nac","momentum")
    include_same_surface_candidates: bool = True
    include_other_surface_candidates: bool = True
    candidate_overlap_block: float = 0.999999

    locality_enter_overlap: float = 0.05
    locality_exit_overlap: float = 0.025

    check_initial_defect: bool = False

    def validate(self):
        if self.defect_interval<=0:
            raise ValueError("defect_interval must be positive.")
        if self.enrich_relative_threshold<=self.prune_relative_threshold:
            raise ValueError("add threshold must exceed prune threshold.")
        if self.min_basis<1 or self.max_basis<self.min_basis:
            raise ValueError("invalid min/max basis.")
        if not (0.0<=self.minimum_capture_fraction<=1.0):
            raise ValueError("minimum_capture_fraction must be in [0,1].")
        if self.minimum_local_utility<0.0:
            raise ValueError("minimum_local_utility cannot be negative.")
        LocalityGraphSettings(
            self.locality_enter_overlap,
            self.locality_exit_overlap,
        ).validate()
        return self


def _copy_basis(basis):
    return [
        DynamicGraphTBF(
            uid=int(b.uid),
            state=int(b.state),
            q=np.asarray(b.q,float).copy(),
            p=np.asarray(b.p,float).copy(),
            A=np.asarray(b.A,float).copy(),
            node=b.node,
            spawned_targets=set(
                getattr(b,"spawned_targets",set())
            ),
        )
        for b in basis
    ]


def _midpoint_basis(old_basis,new_basis):
    if len(old_basis)!=len(new_basis):
        raise ValueError("midpoint basis size changed inside a time step.")
    out=[]
    for old,new in zip(old_basis,new_basis):
        if old.uid!=new.uid:
            raise ValueError("basis uid changed inside a time step.")
        out.append(
            DynamicGraphTBF(
                uid=int(old.uid),
                state=int(old.state),
                q=0.5*(old.q+new.q),
                p=0.5*(old.p+new.p),
                A=np.asarray(old.A,float).copy(),
                node=old.node,
                spawned_targets=set(
                    getattr(old,"spawned_targets",set())
                ),
            )
        )
    return out


def _kinematic_arrays(basis,provider):
    qdots=[]; pdots=[]
    for b in basis:
        qdot,pdot=_kinematics(b,provider)
        qdots.append(qdot); pdots.append(pdot)
    return np.asarray(qdots,float),np.asarray(pdots,float)


def _grid_inner(a,b,area):
    return np.vdot(
        np.asarray(a,dtype=complex).reshape(-1),
        np.asarray(b,dtype=complex).reshape(-1),
    )*float(area)


def _cache_stats(cache):
    return {
        "canonical_solves":int(cache.stats.canonical_solves),
        "requests":int(cache.stats.requests),
    }


def compute_sparse_tdse_defect_v16(
    C,
    basis,
    provider,
    grid,
    mats,
    locality_update,
    ledger=None,
):
    """Evaluate the physical-grid TDSE defect using the sparse projected equation."""
    C=np.asarray(C,dtype=complex)
    qdots,pdots=_kinematic_arrays(basis,provider)

    if ledger is not None:
        ledger.sparse_time_builds+=1

    if ledger is None:
        T=build_sparse_spinor_time_matrix(
            locality_update,qdots,pdots
        )
    else:
        with ledger.timed("time_matrix"):
            T=build_sparse_spinor_time_matrix(
                locality_update,qdots,pdots
            )

    if ledger is not None:
        ledger.sparse_defect_solves+=1

    Cdot=np.asarray(
        spsolve(
            mats.S.tocsc(),
            -(1j*mats.H+T)@C,
        ),
        dtype=complex,
    )

    psi=reconstruct_spinor_complete_wavefunction(
        C,basis,grid.points
    )
    psidot=reconstruct_spinor_complete_time_derivative(
        C,Cdot,basis,grid.points,
        qdots,pdots,
    )
    Hpsi=apply_lvc_grid_hamiltonian(
        psi,grid
    )
    residual=1j*psidot-Hpsi

    r2=float(max(
        np.real(
            _grid_inner(
                residual,residual,grid.area
            )
        ),
        0.0,
    ))
    h2=float(max(
        np.real(
            _grid_inner(
                Hpsi,Hpsi,grid.area
            )
        ),
        0.0,
    ))

    n=len(basis)
    b=np.zeros((n,2),dtype=complex)
    for i,tbf in enumerate(basis):
        g=gaussian_nd(
            grid.points,tbf.q,tbf.p,tbf.A
        )
        for a in range(2):
            b[i,a]=np.vdot(
                g,residual[...,a]
            )*grid.area

    projected=0.0
    Sn=mats.Snuc.tocsc()
    for a in range(2):
        try:
            coeff=np.asarray(
                spsolve(Sn,b[:,a]),
                dtype=complex,
            )
            if not np.all(np.isfinite(coeff)):
                raise FloatingPointError
        except Exception:
            coeff=np.linalg.lstsq(
                mats.Snuc.toarray(),
                b[:,a],
                rcond=1e-12,
            )[0]
        projected+=float(
            np.real(np.vdot(b[:,a],coeff))
        )

    return TDSEDefect(
        residual=residual,
        wavefunction=psi,
        wavefunction_time_derivative=psidot,
        hamiltonian_wavefunction=Hpsi,
        residual_norm=float(np.sqrt(r2)),
        relative_to_hpsi=float(
            np.sqrt(r2)/max(np.sqrt(h2),1e-30)
        ),
        coefficient_derivative=Cdot,
        projected_residual_norm=float(
            np.sqrt(max(projected,0.0))
        ),
    )


def run_sparse_cost_aware_lvc_gaussians(
    initial_basis,
    C0=None,
    provider=None,
    grid=None,
    dt=0.005,
    steps=120,
    settings=SparseAdaptiveSettingsV16(),
    electronic_cost_model=None,
    store_every=10,
):
    """v0.16 sparse/local TDSE-defect-controlled Gaussian propagation."""
    settings=settings.validate()
    provider=provider or AnalyticCI2DFrameProvider()
    if not isinstance(provider,AnalyticCI2DFrameProvider):
        raise TypeError(
            "v0.16 release runner is validated on the analytic 2-state LVC provider."
        )
    if grid is None:
        raise ValueError("diagnostic grid is required.")

    electronic_cost_model=(
        electronic_cost_model
        if electronic_cost_model is not None
        else UniformElectronicCostModel(0.0)
    )

    basis=_copy_basis(initial_basis)
    if not basis:
        raise ValueError("initial basis cannot be empty.")

    if C0 is None:
        C=initialize_spinor_complete_coefficients(
            basis,provider
        )
    else:
        arr=np.asarray(C0,dtype=complex)
        if arr.shape==(len(basis),2):
            C=flatten_coefficients(arr)
        elif arr.shape==(2*len(basis),):
            C=arr.copy()
        elif arr.shape==(len(basis),):
            C=initialize_spinor_complete_coefficients(
                basis,provider,arr
            )
        else:
            raise ValueError("C0 has incompatible shape.")

    graph_settings=LocalityGraphSettings(
        enter_overlap=settings.locality_enter_overlap,
        exit_overlap=settings.locality_exit_overlap,
    )
    endpoint_graph=PersistentGaussianLocalityGraph(
        graph_settings
    )
    midpoint_graph=PersistentGaussianLocalityGraph(
        graph_settings
    )

    ledger=SparseComplexityLedgerV16()
    ledger.start()

    def endpoint_build(cache=None):
        ledger.endpoint_graph_updates+=1
        if cache is None:
            before={"canonical_solves":0,"requests":0}
        else:
            before=_cache_stats(cache)

        with ledger.timed("graph"):
            update=endpoint_graph.update(
                basis,cache=cache
            )
        ledger.observe_graph(update,len(basis))
        ledger.record_pair_delta(
            update.cache,before,
            category="propagation",
        )

        before_matrix=_cache_stats(update.cache)
        ledger.sparse_matrix_builds+=1
        with ledger.timed("matrix"):
            mats=build_sparse_spinor_lvc_matrices(
                update,provider
            )
        ledger.record_pair_delta(
            update.cache,before_matrix,
            category="propagation",
        )
        ledger.observe_matrices(mats)
        return update,mats

    current_update,mats=endpoint_build()

    C=C/np.sqrt(
        sparse_generalized_norm(C,mats.S)
    )

    next_uid=max(int(b.uid) for b in basis)+1
    birth_step={int(b.uid):0 for b in basis}
    lineage={
        int(b.uid):{
            "parent_uid":None,
            "birth_step":0,
            "source":"initial",
            "guidance_state":int(b.state),
        }
        for b in basis
    }

    events=[]
    records=[]
    defect_history=[]
    cost_history=[]
    latest_defect=None
    low_defect_streak=0
    last_adaptation_step=-10**9
    basis_size_sum=0.0

    def current_condition():
        return float(
            np.linalg.cond(
                mats.Snuc.toarray()
            )
        )

    def renormalize():
        nonlocal C
        n=sparse_generalized_norm(C,mats.S)
        if n<=0.0 or not np.isfinite(n):
            raise RuntimeError("invalid sparse generalized norm.")
        C=C/np.sqrt(n)

    def protected_uids(step):
        return {
            uid for uid,born in birth_step.items()
            if step-born<settings.minimum_prune_age_steps
        }

    def evaluate_defect(step):
        ledger.sparse_defect_solves+=0
        with ledger.timed("defect"):
            d=compute_sparse_tdse_defect_v16(
                C,basis,provider,grid,
                mats,current_update,
                ledger=ledger,
            )
        defect_history.append({
            "step":int(step),
            "time":float(step*dt),
            "basis_size":int(len(basis)),
            "relative_to_hpsi":
                float(d.relative_to_hpsi),
            "residual_norm":
                float(d.residual_norm),
            "projected_residual_norm":
                float(d.projected_residual_norm),
            "active_edges":
                int(current_update.active_offdiagonal_edges),
            "edge_fraction":
                float(current_update.edge_fraction),
            "condition_number":
                current_condition(),
        })
        return d

    def pruning_attempt(step,reason,max_loss,require_condition_improvement=False):
        nonlocal basis,C,current_update,mats,last_adaptation_step

        if len(basis)<=settings.min_basis:
            return False

        with ledger.timed("pruning"):
            result=prune_low_loss_gaussian_pair(
                coefficients_matrix(C,len(basis)),
                mats.Snuc.toarray(),
                uids=[b.uid for b in basis],
                max_fractional_loss=max_loss,
                protected_uids=protected_uids(step),
                require_condition_improvement=
                    require_condition_improvement,
            )
        if result is None:
            return False

        removed_uid=int(result.removed_uid)
        keep=result.keep

        basis=[basis[i] for i in keep]
        C=flatten_coefficients(
            result.coefficients_matrix
        )

        subset_cache=current_update.cache.subset(
            keep
        )
        current_update,mats=endpoint_build(
            cache=subset_cache
        )
        renormalize()

        birth_step.pop(removed_uid,None)
        ledger.pruning_events+=1
        last_adaptation_step=int(step)

        events.append({
            "kind":"sparse_prune",
            "step":int(step),
            "time":float(step*dt),
            "reason":str(reason),
            "removed_uid":removed_uid,
            "fractional_projection_loss":
                float(result.fractional_projection_loss),
            "basis_after":int(len(basis)),
            "active_edges_after":
                int(current_update.active_offdiagonal_edges),
        })
        return True

    def enrichment_attempt(step,defect):
        nonlocal basis,C,current_update,mats,next_uid,last_adaptation_step

        dynamic=generate_energy_conserving_defect_candidates_v15(
            basis,
            provider,
            position_shifts=
                settings.candidate_position_shifts,
            width_scales=
                settings.candidate_width_scales,
            momentum_directions=
                settings.candidate_momentum_directions,
            include_same_surface=
                settings.include_same_surface_candidates,
            include_other_surfaces=
                settings.include_other_surface_candidates,
            overlap_block=
                settings.candidate_overlap_block,
        )

        ledger.candidate_searches+=1
        ledger.candidates_scored+=len(dynamic)

        with ledger.timed("candidate"):
            ranked=rank_dynamic_defect_candidates_cached(
                defect,
                basis,
                dynamic,
                grid,
                current_update.cache,
                mats.Snuc.toarray(),
                condition_limit=settings.condition_limit,
                orthogonal_norm_floor=
                    settings.orthogonal_norm_floor,
                exact_condition_top=
                    settings.residual_shortlist,
                max_return=
                    settings.residual_shortlist,
            )

        for score in ranked:
            temp=score.expanded_cache
            ledger.pair_factorizations+=int(
                temp.stats.canonical_solves
            )
            ledger.candidate_pair_factorizations+=int(
                temp.stats.canonical_solves
            )
            ledger.pair_requests+=int(
                temp.stats.requests
            )

        if not ranked:
            events.append({
                "kind":"sparse_enrichment_rejected",
                "step":int(step),
                "reason":"no_residual_candidate",
            })
            return False,None

        ledger.cost_reranks+=1
        with ledger.timed("cost"):
            cost_ranked=rank_local_sparse_candidates(
                ranked,
                dynamic,
                basis,
                active_offdiagonal_edges=
                    current_update.active_offdiagonal_edges,
                overlap_threshold=
                    settings.locality_enter_overlap,
                current_condition=
                    current_condition(),
                horizon_steps=
                    settings.cost_horizon_steps,
                minimum_capture_fraction=
                    settings.minimum_capture_fraction,
                minimum_utility=
                    settings.minimum_local_utility,
                condition_penalty_weight=
                    settings.condition_penalty_weight,
                electronic_cost_model=
                    electronic_cost_model,
                electronic_cost_weight=
                    settings.electronic_cost_weight,
            )

        if not cost_ranked:
            events.append({
                "kind":"sparse_enrichment_rejected",
                "step":int(step),
                "reason":"local_cost_gate",
                "residual_best_capture":
                    float(ranked[0].capture_fraction),
            })
            return False,None

        best_cost=cost_ranked[0]
        best_by_index={
            x.candidate_index:x
            for x in ranked
        }
        best=best_by_index[
            best_cost.candidate_index
        ]
        item=dynamic[
            best_cost.candidate_index
        ]

        child=item.candidate.to_tbf(
            next_uid,
            node_prefix="v16_sparse",
        )
        estimate=electronic_cost_model.estimate(
            child.q
        )
        ledger.record_electronic_cost(
            estimate
        )

        basis_before=len(basis)
        basis.append(child)
        C=np.concatenate([
            C,
            np.zeros(2,dtype=complex),
        ])

        expanded_cache=best.expanded_cache
        expanded_cache.basis[-1]=child

        current_update,mats=endpoint_build(
            cache=expanded_cache
        )
        renormalize()

        electronic_cost_model.register(
            child.q
        )

        birth_step[int(next_uid)]=int(step)
        lineage[int(next_uid)]={
            "parent_uid":int(item.parent_uid),
            "birth_step":int(step),
            "source":"sparse_cost_aware_tdse_defect",
            "guidance_state":int(child.state),
            "candidate_label":
                str(item.candidate.label),
        }
        next_uid+=1
        last_adaptation_step=int(step)
        ledger.enrichments+=1

        after=evaluate_defect(step)

        cost_history.append({
            "step":int(step),
            "candidate_label":
                str(item.candidate.label),
            "capture_fraction":
                float(best.capture_fraction),
            "utility":float(best_cost.utility),
            "normalized_incremental_cost":
                float(best_cost.normalized_incremental_cost),
            "predicted_local_degree":
                int(best_cost.local_degree),
            "estimated_nnz_growth":
                int(best_cost.estimated_nnz_growth),
            "electronic_cost_units":
                float(best_cost.electronic_cost_units),
            "electronic_cache_hit":
                bool(best_cost.electronic_cache_hit),
        })

        events.append({
            "kind":"sparse_cost_aware_enrichment",
            "step":int(step),
            "time":float(step*dt),
            "parent_uid":int(item.parent_uid),
            "new_uid":int(child.uid),
            "candidate_label":
                str(item.candidate.label),
            "basis_before":int(basis_before),
            "basis_after":int(len(basis)),
            "relative_defect_before":
                float(defect.relative_to_hpsi),
            "relative_defect_after":
                float(after.relative_to_hpsi),
            "capture_fraction_predicted":
                float(best.capture_fraction),
            "utility":float(best_cost.utility),
            "normalized_incremental_cost":
                float(best_cost.normalized_incremental_cost),
            "predicted_local_degree":
                int(best_cost.local_degree),
            "electronic_cost_units":
                float(best_cost.electronic_cost_units),
            "electronic_cache_hit":
                bool(best_cost.electronic_cache_hit),
            "active_edges_after":
                int(current_update.active_offdiagonal_edges),
            "edge_fraction_after":
                float(current_update.edge_fraction),
            "zero_coefficient_insertion":True,
        })
        return True,after

    def control(step):
        nonlocal latest_defect,low_defect_streak,last_adaptation_step

        latest_defect=evaluate_defect(step)

        if (
            current_condition()>settings.hard_condition_limit
            and len(basis)>settings.min_basis
        ):
            if pruning_attempt(
                step,
                "hard_condition_limit",
                settings.emergency_prune_fractional_loss,
                True,
            ):
                latest_defect=evaluate_defect(step)

        separated=(
            step-last_adaptation_step
            >=settings.minimum_adaptation_separation_steps
        )

        if latest_defect.relative_to_hpsi>=settings.enrich_relative_threshold:
            low_defect_streak=0
            if not separated:
                return

            if len(basis)>=settings.max_basis:
                if pruning_attempt(
                    step,
                    "basis_budget_replacement",
                    settings.max_replacement_prune_fractional_loss,
                    False,
                ):
                    latest_defect=evaluate_defect(step)

            if len(basis)<settings.max_basis:
                changed,after=enrichment_attempt(
                    step,latest_defect
                )
                if changed:
                    latest_defect=after

        elif latest_defect.relative_to_hpsi<=settings.prune_relative_threshold:
            low_defect_streak+=1
            if (
                low_defect_streak>=settings.prune_patience_checks
                and separated
                and len(basis)>settings.min_basis
            ):
                if pruning_attempt(
                    step,
                    "sustained_low_defect",
                    settings.max_prune_fractional_loss,
                    False,
                ):
                    low_defect_streak=0
                    latest_defect=evaluate_defect(step)
        else:
            low_defect_streak=0

    def record(step):
        rho=sparse_reduced_density(
            C,mats.Snuc,normalize=True
        )
        records.append({
            "step":int(step),
            "time":float(step*dt),
            "norm":float(
                sparse_generalized_norm(
                    C,mats.S
                )
            ),
            "basis_size":int(len(basis)),
            "active_edges":
                int(current_update.active_offdiagonal_edges),
            "edge_fraction":
                float(current_update.edge_fraction),
            "S_nnz":int(mats.S.nnz),
            "H_nnz":int(mats.H.nnz),
            "condition_number":
                current_condition(),
            "diabatic_populations":
                np.real(np.diag(rho)),
            "diabatic_coherence":
                complex(rho[0,1]),
            "latest_relative_defect":
                None if latest_defect is None
                else float(latest_defect.relative_to_hpsi),
        })

    if settings.check_initial_defect:
        control(0)
    record(0)

    for step in range(1,int(steps)+1):
        old_basis=_copy_basis(basis)
        old_mats=mats

        qdot0,pdot0=_kinematic_arrays(
            old_basis,provider
        )
        for b in basis:
            b.q,b.p=_verlet(
                b,provider,dt
            )
        qdot1,pdot1=_kinematic_arrays(
            basis,provider
        )

        current_update,mats=endpoint_build()

        mid=_midpoint_basis(
            old_basis,basis
        )
        qdot_mid=0.5*(qdot0+qdot1)
        pdot_mid=0.5*(pdot0+pdot1)

        ledger.midpoint_graph_updates+=1
        with ledger.timed("graph"):
            mid_update=midpoint_graph.update(mid)
        ledger.observe_graph(
            mid_update,len(mid)
        )
        # Fresh midpoint cache: graph exact checks are propagation pair work.
        ledger.record_pair_delta(
            mid_update.cache,
            {"canonical_solves":0,"requests":0},
            category="propagation",
        )

        before_mid=_cache_stats(mid_update.cache)
        ledger.sparse_time_builds+=1
        with ledger.timed("time_matrix"):
            Tseed=build_sparse_spinor_time_matrix(
                mid_update,
                qdot_mid,
                pdot_mid,
            )
        ledger.record_pair_delta(
            mid_update.cache,
            before_mid,
            category="propagation",
        )

        Tmid=sparse_metric_compatible_connection(
            old_mats.S,
            mats.S,
            dt,
            Tseed,
        )

        ledger.sparse_cayley_solves+=1
        with ledger.timed("cayley"):
            C=sparse_moving_basis_midpoint_cayley_step(
                C,
                old_mats.S,
                old_mats.H,
                mats.S,
                mats.H,
                Tmid,
                dt,
            )

        if step%settings.defect_interval==0:
            control(step)

        ledger.peak_basis_size=max(
            ledger.peak_basis_size,
            len(basis),
        )
        basis_size_sum+=len(basis)

        if step%int(store_every)==0:
            record(step)

    if records[-1]["step"]!=int(steps):
        record(int(steps))

    ledger.stop()

    return {
        "records":records,
        "defect_history":defect_history,
        "cost_history":cost_history,
        "events":events,
        "final_basis":basis,
        "final_coefficients":C,
        "final_sparse_matrices":mats,
        "lineage":lineage,
        "average_basis_size":
            float(basis_size_sum/max(int(steps),1)),
        "endpoint_graph":
            endpoint_graph.diagnostics(),
        "midpoint_graph":
            midpoint_graph.diagnostics(),
        "complexity":ledger.as_dict(),
        "settings":{
            "dt":float(dt),
            "steps":int(steps),
            "store_every":int(store_every),
            "control":asdict(settings),
            "representation":
                "local-overlap sparse spinor-complete Gaussian basis",
            "solver":"scipy.sparse.linalg.spsolve",
            "locality":
                "conservative overlap-bound screen + exact overlap + hysteresis",
            "adaptation":
                "TDSE residual benefit / local sparse + electronic cost",
        },
    }
