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
)
from .defect_candidates_v18 import (
    rank_dynamic_defect_candidates_batched_v18,
)
from .edge_importance_v17 import (
    EdgeImportanceSettingsV17,
    ErrorControlledGaussianLocalityGraphV17,
)
from .sparse_pair_matrices_v16 import (
    build_sparse_spinor_lvc_matrices,
    build_sparse_spinor_time_matrix,
    sparse_metric_compatible_connection,
    sparse_moving_basis_midpoint_cayley_step,
    sparse_generalized_norm,
    sparse_reduced_density,
    audit_sparse_lvc_matrices_against_dense,
)
from .local_cost_aware_v16 import (
    rank_local_sparse_candidates,
)
from .electronic_cost_v16 import (
    UniformElectronicCostModel,
    ElectronicCostEstimate,
)
from .sampled_sparse_audit_v18 import (
    sampled_omitted_edge_audit_v18,
)
from .convergence_complexity_v18 import (
    ConvergenceComplexityLedgerV18,
)


@dataclass(frozen=True)
class ConvergenceCompleteSettingsV18:
    defect_interval: int = 10
    defect_interval_time: float | None = None
    enrich_relative_threshold: float = 0.020
    prune_relative_threshold: float = 0.006

    minimum_capture_fraction: float = 0.003
    minimum_local_utility: float = 0.08
    condition_penalty_weight: float = 0.15
    electronic_cost_weight: float = 1.0
    cost_horizon_steps: int = 10
    cost_horizon_time: float | None = None
    residual_shortlist: int = 8

    min_basis: int = 8
    max_basis: int = 11
    minimum_adaptation_separation_steps: int = 10
    minimum_adaptation_separation_time: float | None = None
    minimum_prune_age_steps: int = 20
    minimum_prune_age_time: float | None = None
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

    edge_enter_score: float = 0.030
    edge_exit_score: float = 0.015
    search_overlap_floor: float = 1e-4
    edge_overlap_weight: float = 1.0
    edge_hamiltonian_weight: float = 0.20
    edge_time_connection_weight: float = 1.0
    local_omitted_score_l2_budget: float = 0.08

    sampled_audit_interval: int = 20
    sampled_audit_interval_time: float | None = None
    sampled_audit_priority_pairs: int = 8
    sampled_audit_random_pairs: int = 8
    sampled_audit_wider_search_factor: float = 0.1
    sampled_audit_seed: int = 20260813
    sampled_audit_violation_factor: float = 1.0
    sampled_audit_relaxation_factor: float = 0.5
    max_sampled_audit_relaxations: int = 3

    sentinel_max_S_error: float = 0.006
    sentinel_max_H_error: float = 0.006
    sentinel_max_Snuc_error: float = 0.006

    candidate_batch_size: int = 16

    check_initial_defect: bool = False

    def validate(self):
        if self.defect_interval<=0:
            raise ValueError("defect_interval must be positive.")
        for name,value in (
            ("defect_interval_time",self.defect_interval_time),
            ("minimum_adaptation_separation_time",self.minimum_adaptation_separation_time),
            ("minimum_prune_age_time",self.minimum_prune_age_time),
            ("sampled_audit_interval_time",self.sampled_audit_interval_time),
            ("cost_horizon_time",self.cost_horizon_time),
        ):
            if value is not None and float(value)<=0.0:
                raise ValueError(f"{name} must be positive when provided.")
        if self.enrich_relative_threshold<=self.prune_relative_threshold:
            raise ValueError("add threshold must exceed prune threshold.")
        if self.min_basis<1 or self.max_basis<self.min_basis:
            raise ValueError("invalid min/max basis.")
        if not (0.0<=self.minimum_capture_fraction<=1.0):
            raise ValueError("minimum_capture_fraction must be in [0,1].")
        if self.minimum_local_utility<0.0:
            raise ValueError("minimum_local_utility cannot be negative.")
        EdgeImportanceSettingsV17(
            enter_score=self.edge_enter_score,
            exit_score=self.edge_exit_score,
            search_overlap_floor=self.search_overlap_floor,
            overlap_weight=self.edge_overlap_weight,
            hamiltonian_weight=self.edge_hamiltonian_weight,
            time_connection_weight=self.edge_time_connection_weight,
            local_omitted_score_l2_budget=
                self.local_omitted_score_l2_budget,
        ).validate()
        if self.sampled_audit_interval<=0:
            raise ValueError("sampled_audit_interval must be positive.")
        if self.sampled_audit_priority_pairs<0 or self.sampled_audit_random_pairs<0:
            raise ValueError("sampled audit counts cannot be negative.")
        if not (0.0<self.sampled_audit_wider_search_factor<=1.0):
            raise ValueError("sampled_audit_wider_search_factor must lie in (0,1].")
        if self.sampled_audit_violation_factor<=0.0:
            raise ValueError("sampled_audit_violation_factor must be positive.")
        if not (0.0<self.sampled_audit_relaxation_factor<1.0):
            raise ValueError("sampled_audit_relaxation_factor must lie in (0,1).")
        if self.max_sampled_audit_relaxations<0:
            raise ValueError("max_sampled_audit_relaxations cannot be negative.")
        if self.candidate_batch_size<=0:
            raise ValueError("candidate_batch_size must be positive.")
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


def run_convergence_complete_lvc_gaussians(
    initial_basis,
    C0=None,
    provider=None,
    grid=None,
    dt=0.005,
    steps=120,
    settings=ConvergenceCompleteSettingsV18(),
    electronic_cost_model=None,
    store_every=10,
    return_snapshots=False,
):
    """v0.18 sampled-audit, batched-candidate sparse Gaussian propagation."""
    settings=settings.validate()
    provider=provider or AnalyticCI2DFrameProvider()
    if not isinstance(provider,AnalyticCI2DFrameProvider):
        raise TypeError(
            "v0.18 release runner is validated on the analytic 2-state LVC provider."
        )
    if grid is None:
        raise ValueError("diagnostic grid is required.")

    electronic_cost_model=(
        electronic_cost_model
        if electronic_cost_model is not None
        else UniformElectronicCostModel(0.0)
    )

    def _resolved_steps(time_value,step_value):
        if time_value is None:
            return max(int(step_value),1)
        return max(int(round(float(time_value)/float(dt))),1)

    defect_interval_steps=_resolved_steps(
        settings.defect_interval_time,
        settings.defect_interval,
    )
    sampled_audit_interval_steps=_resolved_steps(
        settings.sampled_audit_interval_time,
        settings.sampled_audit_interval,
    )
    adaptation_separation_steps=_resolved_steps(
        settings.minimum_adaptation_separation_time,
        settings.minimum_adaptation_separation_steps,
    )
    prune_age_steps=_resolved_steps(
        settings.minimum_prune_age_time,
        settings.minimum_prune_age_steps,
    )
    cost_horizon_steps=_resolved_steps(
        settings.cost_horizon_time,
        settings.cost_horizon_steps,
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

    graph_settings=EdgeImportanceSettingsV17(
        enter_score=settings.edge_enter_score,
        exit_score=settings.edge_exit_score,
        search_overlap_floor=settings.search_overlap_floor,
        overlap_weight=settings.edge_overlap_weight,
        hamiltonian_weight=settings.edge_hamiltonian_weight,
        time_connection_weight=settings.edge_time_connection_weight,
        local_omitted_score_l2_budget=
            settings.local_omitted_score_l2_budget,
    )
    endpoint_graph=ErrorControlledGaussianLocalityGraphV17(
        provider,dt,graph_settings
    )
    midpoint_graph=ErrorControlledGaussianLocalityGraphV17(
        provider,dt,graph_settings
    )

    ledger=ConvergenceComplexityLedgerV18()
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
    sampled_audit_history=[]
    sentinel_audit_history=[]
    candidate_batch_history=[]
    snapshots=[]
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
            if step-born<prune_age_steps
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
            ranked,batch_diag=rank_dynamic_defect_candidates_batched_v18(
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
                batch_size=settings.candidate_batch_size,
                return_diagnostics=True,
            )
        ledger.record_candidate_batching(batch_diag)
        candidate_batch_history.append({
            "step":int(step),
            "time":float(step*dt),
            **batch_diag.as_dict(),
        })

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
                    settings.search_overlap_floor,
                current_condition=
                    current_condition(),
                horizon_steps=
                    cost_horizon_steps,
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
            >=adaptation_separation_steps
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

    def sentinel_dense_audit(step,label):
        """Full dense S/H/Snuc audit reserved for initial/final release sentinels."""
        with ledger.timed("audit"):
            audit=audit_sparse_lvc_matrices_against_dense(
                basis,provider,mats
            )
        ledger.record_sentinel_audit(audit)

        passed=(
            audit["relative_S_frobenius_error"]
            <=settings.sentinel_max_S_error
            and audit["relative_H_frobenius_error"]
            <=settings.sentinel_max_H_error
            and audit["relative_Snuc_frobenius_error"]
            <=settings.sentinel_max_Snuc_error
        )
        row={
            "step":int(step),
            "time":float(step*dt),
            "label":str(label),
            "passed":bool(passed),
            "enter_score":
                float(endpoint_graph.settings.enter_score),
            "exit_score":
                float(endpoint_graph.settings.exit_score),
            "search_overlap_floor":
                float(endpoint_graph.settings.search_overlap_floor),
            **audit,
        }
        sentinel_audit_history.append(row)
        if not passed:
            events.append({
                "kind":"sentinel_dense_audit_failed",
                **row,
            })
        return row

    def sampled_audit_and_relax(step):
        """Normal v0.18 audit: exact S/H/T scoring on a deterministic omitted-edge sample."""
        nonlocal current_update,mats

        if step%sampled_audit_interval_steps!=0:
            return None

        attempts=0
        while True:
            with ledger.timed("sampled_audit"):
                audit=sampled_omitted_edge_audit_v18(
                    basis,
                    provider,
                    dt,
                    current_update,
                    endpoint_graph.settings,
                    step=step,
                    priority_count=
                        settings.sampled_audit_priority_pairs,
                    random_count=
                        settings.sampled_audit_random_pairs,
                    wider_search_factor=
                        settings.sampled_audit_wider_search_factor,
                    seed=settings.sampled_audit_seed,
                    violation_factor=
                        settings.sampled_audit_violation_factor,
                )
            ledger.record_sampled_audit(audit)

            row={
                **audit.as_dict(),
                "attempt":int(attempts),
                "enter_score":
                    float(endpoint_graph.settings.enter_score),
                "exit_score":
                    float(endpoint_graph.settings.exit_score),
                "search_overlap_floor":
                    float(endpoint_graph.settings.search_overlap_floor),
            }
            sampled_audit_history.append(row)

            if audit.passed:
                return row

            if attempts>=settings.max_sampled_audit_relaxations:
                events.append({
                    "kind":"sampled_sparse_audit_unresolved",
                    **row,
                })
                return row

            factor=settings.sampled_audit_relaxation_factor
            endpoint_graph.relax_scores(factor)
            midpoint_graph.relax_scores(factor)
            endpoint_graph.relax_search_floor(factor)
            midpoint_graph.relax_search_floor(factor)
            ledger.score_relaxations+=1
            ledger.search_floor_relaxations+=1

            events.append({
                "kind":"sampled_sparse_audit_relaxation",
                "step":int(step),
                "time":float(step*dt),
                "attempt":int(attempts),
                "maximum_sampled_score":
                    float(audit.maximum_score),
                "sampled_violation_count":
                    int(audit.violation_count),
                "new_enter_score":
                    float(endpoint_graph.settings.enter_score),
                "new_exit_score":
                    float(endpoint_graph.settings.exit_score),
                "new_search_overlap_floor":
                    float(endpoint_graph.settings.search_overlap_floor),
            })

            current_update,mats=endpoint_build()
            renormalize()
            attempts+=1

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
            "omitted_candidate_score_l2":
                float(current_update.omitted_candidate_score_l2),
            "omitted_candidate_score_max":
                float(current_update.omitted_candidate_score_max),
            "budget_promoted_edges":
                int(current_update.budget_promoted_edges),
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
        if return_snapshots:
            snapshots.append({
                "step":int(step),
                "time":float(step*dt),
                "basis":_copy_basis(basis),
                "coefficients":np.asarray(C,dtype=complex).copy(),
            })

    initial_sentinel=sentinel_dense_audit(0,"initial")
    if not initial_sentinel["passed"]:
        raise RuntimeError(
            "v0.18 initial sparse sentinel failed; calibrate graph settings before propagation."
        )
    renormalize()

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
        sampled_audit_and_relax(step)

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

        if step%defect_interval_steps==0:
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

    final_sentinel=sentinel_dense_audit(
        int(steps),"final"
    )
    renormalize()

    ledger.stop()

    return {
        "records":records,
        "defect_history":defect_history,
        "cost_history":cost_history,
        "sampled_audit_history":
            sampled_audit_history,
        "sentinel_audit_history":
            sentinel_audit_history,
        "candidate_batch_history":
            candidate_batch_history,
        "snapshots":snapshots,
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
            "resolved_control_steps":{
                "defect_interval":int(defect_interval_steps),
                "sampled_audit_interval":int(sampled_audit_interval_steps),
                "adaptation_separation":int(adaptation_separation_steps),
                "prune_age":int(prune_age_steps),
                "cost_horizon":int(cost_horizon_steps),
            },
            "control":asdict(settings),
            "representation":
                "sampled-audited S/H/T error-controlled sparse spinor-complete Gaussian basis",
            "solver":"scipy.sparse.linalg.spsolve",
            "locality":
                "safe geometric pre-screen + exact local S/H/T score + hysteresis + sampled omitted-edge audits + initial/final dense sentinels",
            "adaptation":
                "TDSE residual benefit / local sparse + electronic cost with batched KxG contractions",
        },
    }
