from dataclasses import dataclass, asdict
import math
import numpy as np

from .dynamic_gauge_graph import AnalyticCI2DFrameProvider
from .dynamic_graph_aims import DynamicGraphTBF, _verlet, _kinematics
from .moving_graph_gaussian import metric_compatible_basis_connection
from .moving_basis_v12 import moving_basis_midpoint_cayley_step
from .spinor_complete_lvc_v12 import (
    spinor_complete_generalized_norm,
    coefficients_matrix,
    flatten_coefficients,
)
from .spinor_complete_dynamics_v12 import (
    initialize_spinor_complete_coefficients,
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
from .pair_cache_v15 import (
    GaussianPairCache,
    build_cached_spinor_lvc_matrices,
    build_cached_spinor_time_matrix,
    expand_cached_spinor_lvc_matrices,
    subset_cached_spinor_lvc_matrices,
    v14_factorization_equivalent_for_sh,
    v14_factorization_equivalent_for_time,
)
from .defect_candidates_v15 import (
    generate_energy_conserving_defect_candidates_v15,
    rank_dynamic_defect_candidates_cached,
)
from .cost_aware_adaptation_v15 import (
    rank_candidates_by_cost_aware_utility,
)
from .complexity_v15 import (
    ComplexityLedgerV15,
    dense_cubic_units,
)


@dataclass(frozen=True)
class AdaptiveDefectSettingsV15:
    """v0.15 residual + computational-cost control policy."""

    defect_interval: int = 10
    enrich_relative_threshold: float = 0.020
    prune_relative_threshold: float = 0.006

    minimum_capture_fraction: float = 0.003
    minimum_cost_aware_utility: float = 0.0
    condition_penalty_weight: float = 0.15
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

    check_initial_defect: bool = False

    def validate(self):
        if self.defect_interval<=0:
            raise ValueError("defect_interval must be positive.")
        if self.enrich_relative_threshold<=self.prune_relative_threshold:
            raise ValueError(
                "enrichment threshold must exceed pruning threshold."
            )
        if self.cost_horizon_steps<=0:
            raise ValueError("cost_horizon_steps must be positive.")
        if self.residual_shortlist<=0:
            raise ValueError("residual_shortlist must be positive.")
        if self.min_basis<1 or self.max_basis<self.min_basis:
            raise ValueError("invalid min/max basis.")
        if not (0.0<=self.minimum_capture_fraction<=1.0):
            raise ValueError("minimum_capture_fraction must be in [0,1].")
        if self.minimum_cost_aware_utility<0.0:
            raise ValueError("minimum_cost_aware_utility cannot be negative.")
        if self.condition_penalty_weight<0.0:
            raise ValueError("condition_penalty_weight cannot be negative.")
        if self.condition_limit<=1.0:
            raise ValueError("condition_limit must exceed one.")
        if self.hard_condition_limit<self.condition_limit:
            raise ValueError(
                "hard_condition_limit must not be smaller than condition_limit."
            )
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
        raise ValueError(
            "midpoint basis requires equal basis sizes."
        )
    out=[]
    for old,new in zip(old_basis,new_basis):
        if old.uid!=new.uid:
            raise ValueError(
                "basis identity changed inside one propagation step."
            )
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
    qdots=[]
    pdots=[]
    for b in basis:
        qdot,pdot=_kinematics(b,provider)
        qdots.append(qdot)
        pdots.append(pdot)
    return np.asarray(qdots,float),np.asarray(pdots,float)


def _grid_inner(a,b,area):
    return np.vdot(
        np.asarray(a,dtype=complex).reshape(-1),
        np.asarray(b,dtype=complex).reshape(-1),
    )*float(area)


def _cache_stat_snapshot(cache):
    s=cache.stats
    return {
        "requests":int(s.requests),
        "canonical_solves":int(s.canonical_solves),
        "direct_hits":int(s.direct_hits),
        "reverse_views":int(s.reverse_views),
        "inherited_pairs":int(s.inherited_pairs),
    }


def reduced_density_from_snuc(
    flat_coefficients,
    Snuc,
    normalize=True,
):
    C=coefficients_matrix(
        flat_coefficients,
        len(Snuc),
    )
    S=np.asarray(Snuc,dtype=complex)

    rho=C.T@S.T@np.conj(C)
    rho=0.5*(rho+rho.conj().T)

    if normalize:
        tr=np.trace(rho)
        if abs(tr)<1e-15:
            raise ValueError("zero reduced-density trace.")
        rho=rho/tr
    return rho


def compute_tdse_defect_cached_v15(
    C,
    basis,
    provider,
    grid,
    S,
    H,
    Snuc,
    cache,
    ledger=None,
):
    """TDSE defect reusing the current endpoint pair cache."""
    C=np.asarray(C,dtype=complex)
    qdots,pdots=_kinematic_arrays(
        basis,provider
    )

    before=_cache_stat_snapshot(cache)
    if ledger is not None:
        ledger.time_matrix_calls+=1
        ledger.v14_factorization_baseline += (
            v14_factorization_equivalent_for_time(
                len(basis)
            )
        )

    if ledger is None:
        T=build_cached_spinor_time_matrix(
            cache,qdots,pdots
        )
    else:
        with ledger.timed("time_matrix"):
            T=build_cached_spinor_time_matrix(
                cache,qdots,pdots
            )
        ledger.add_cache_delta(
            cache.stats,before,
            category="propagation",
        )

    m=len(C)
    if ledger is not None:
        ledger.defect_solve_calls+=1
        ledger.defect_cubic_units+=m**3

    Cdot=np.linalg.solve(
        np.asarray(S,dtype=complex),
        -(1j*np.asarray(H,dtype=complex)+T)@C,
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
    for a in range(2):
        coeff=np.linalg.lstsq(
            Snuc,b[:,a],rcond=1e-12
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


def run_time_adaptive_cost_aware_lvc_gaussians(
    initial_basis,
    C0=None,
    provider=None,
    grid=None,
    dt=0.005,
    steps=120,
    settings=AdaptiveDefectSettingsV15(),
    store_every=10,
):
    r"""v0.15 time-adaptive runner with cached pair algebra and cost-aware growth."""
    settings=settings.validate()
    provider=provider or AnalyticCI2DFrameProvider()

    if not isinstance(
        provider,AnalyticCI2DFrameProvider
    ):
        raise TypeError(
            "v0.15 reference runner is specialized to the analytic 2-state LVC "
            "provider."
        )
    if grid is None:
        raise ValueError(
            "A diagnostic grid is required."
        )
    if int(store_every)<=0:
        raise ValueError(
            "store_every must be positive."
        )

    basis=_copy_basis(initial_basis)
    if not basis:
        raise ValueError(
            "initial basis cannot be empty."
        )

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
            raise ValueError(
                "C0 has incompatible shape."
            )

    ledger=ComplexityLedgerV15()
    ledger.start()

    def full_endpoint_build():
        ledger.full_matrix_builds+=1
        ledger.pair_snapshots+=1
        ledger.v14_factorization_baseline += (
            v14_factorization_equivalent_for_sh(
                len(basis)
            )
        )

        cache=GaussianPairCache(basis)
        before=_cache_stat_snapshot(cache)
        with ledger.timed("matrix_build"):
            S,H,Snuc=build_cached_spinor_lvc_matrices(
                cache,provider
            )
        ledger.add_cache_delta(
            cache.stats,before,
            category="propagation",
        )
        return S,H,Snuc,cache

    S,H,Snuc,cache=full_endpoint_build()
    C/=np.sqrt(
        spinor_complete_generalized_norm(C,S)
    )

    next_uid=max(int(b.uid) for b in basis)+1
    birth_step={
        int(b.uid):0
        for b in basis
    }
    lineage={
        int(b.uid):{
            "parent_uid":None,
            "birth_step":0,
            "birth_time":0.0,
            "guidance_state":int(b.state),
            "source":"initial",
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
    basis_size_samples=0

    latest_connection_correction=0.0
    latest_connection_seed=0.0

    def current_condition():
        return float(np.linalg.cond(Snuc))

    def protected_uids(step):
        return {
            int(uid)
            for uid,born in birth_step.items()
            if (
                step-int(born)
                <settings.minimum_prune_age_steps
            )
        }

    def renormalize():
        nonlocal C
        norm=spinor_complete_generalized_norm(
            C,S
        )
        if norm<=0.0 or not np.isfinite(norm):
            raise RuntimeError(
                "invalid generalized norm."
            )
        C=C/np.sqrt(norm)

    def pruning_attempt(
        step,
        reason,
        max_loss,
        require_condition_improvement=False,
    ):
        nonlocal basis,C,S,H,Snuc,cache
        nonlocal last_adaptation_step

        if len(basis)<=settings.min_basis:
            return False

        ledger.pruning_audits+=1
        with ledger.timed("pruning"):
            Cmat=coefficients_matrix(
                C,len(basis)
            )
            result=prune_low_loss_gaussian_pair(
                Cmat,
                Snuc,
                uids=[b.uid for b in basis],
                max_fractional_loss=max_loss,
                protected_uids=protected_uids(step),
                require_condition_improvement=
                    require_condition_improvement,
            )

        if result is None:
            return False

        old_basis_size=len(basis)
        removed_uid=int(result.removed_uid)

        # v0.14 rebuilt all S/H pair algebra after pruning.
        nnew=len(result.keep)
        ledger.v14_factorization_baseline += (
            v14_factorization_equivalent_for_sh(
                nnew
            )
        )

        S,H,Snuc,cache=(
            subset_cached_spinor_lvc_matrices(
                S,H,Snuc,cache,result.keep
            )
        )
        ledger.incremental_prunes+=1
        ledger.inherited_pairs_reused += (
            cache.stats.inherited_pairs
        )

        basis=[basis[i] for i in result.keep]
        C=flatten_coefficients(
            result.coefficients_matrix
        )
        renormalize()

        birth_step.pop(removed_uid,None)
        lineage.setdefault(
            removed_uid,{}
        )["removed_step"]=int(step)
        lineage[removed_uid][
            "removed_time"
        ]=float(step*dt)
        lineage[removed_uid][
            "removal_reason"
        ]=str(reason)

        ledger.pruning_events+=1
        last_adaptation_step=int(step)

        events.append({
            "kind":"residual_prune_v15",
            "step":int(step),
            "time":float(step*dt),
            "reason":str(reason),
            "removed_uid":removed_uid,
            "basis_before":int(old_basis_size),
            "basis_after":int(len(basis)),
            "absolute_projection_loss":
                float(result.absolute_projection_loss),
            "fractional_projection_loss":
                float(result.fractional_projection_loss),
            "condition_before":
                float(result.condition_before),
            "condition_after":
                float(result.condition_after),
            "incremental_matrix_slice":True,
            "pair_factorizations_for_prune":0,
        })
        return True

    def evaluate_defect(step):
        ledger.defect_evaluations+=1
        with ledger.timed("defect"):
            d=compute_tdse_defect_cached_v15(
                C,basis,provider,grid,
                S,H,Snuc,cache,
                ledger=ledger,
            )

        defect_history.append({
            "step":int(step),
            "time":float(step*dt),
            "basis_size":int(len(basis)),
            "residual_norm":
                float(d.residual_norm),
            "relative_to_hpsi":
                float(d.relative_to_hpsi),
            "projected_residual_norm":
                float(d.projected_residual_norm),
            "condition_number":
                current_condition(),
        })
        return d

    def enrichment_attempt(step,defect):
        nonlocal basis,C,S,H,Snuc,cache
        nonlocal next_uid,last_adaptation_step

        dynamic=(
            generate_energy_conserving_defect_candidates_v15(
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
        )

        ledger.candidate_ranking_calls+=1
        ledger.candidate_count_scored+=len(dynamic)
        ledger.peak_candidate_count=max(
            ledger.peak_candidate_count,
            len(dynamic),
        )

        # Track pair work created by exact condition checks on the residual shortlist.
        candidate_caches_before=[]

        with ledger.timed("candidate_ranking"):
            ranked=rank_dynamic_defect_candidates_cached(
                defect,
                basis,
                dynamic,
                grid,
                cache,
                Snuc,
                condition_limit=settings.condition_limit,
                orthogonal_norm_floor=
                    settings.orthogonal_norm_floor,
                exact_condition_top=
                    settings.residual_shortlist,
                max_return=
                    settings.residual_shortlist,
            )

        # Each returned score owns a temporary expanded cache.  Its old-old pairs are
        # inherited and only child pairs were factorized.
        for score in ranked:
            temp=score.expanded_cache
            ledger.pair_factorizations += (
                temp.stats.canonical_solves
            )
            ledger.candidate_pair_factorizations += (
                temp.stats.canonical_solves
            )
            ledger.pair_requests += (
                temp.stats.requests
            )
            ledger.pair_direct_hits += (
                temp.stats.direct_hits
            )
            ledger.pair_reverse_views += (
                temp.stats.reverse_views
            )
            ledger.inherited_pairs_reused += (
                temp.stats.inherited_pairs
            )

        if not ranked:
            events.append({
                "kind":"enrichment_rejected_v15",
                "step":int(step),
                "time":float(step*dt),
                "reason":"no_admissible_candidate",
                "basis_size":int(len(basis)),
                "relative_defect":
                    float(defect.relative_to_hpsi),
                "candidate_count":int(len(dynamic)),
            })
            return False,None

        rates=ledger.empirical_cost_rates()
        defect_checks=max(
            1,
            math.ceil(
                settings.cost_horizon_steps
                /settings.defect_interval
            ),
        )

        ledger.cost_ranking_calls+=1
        with ledger.timed("cost_ranking"):
            cost_ranked=rank_candidates_by_cost_aware_utility(
                ranked,
                n_basis=len(basis),
                current_condition=current_condition(),
                horizon_steps=
                    settings.cost_horizon_steps,
                defect_checks=defect_checks,
                minimum_capture_fraction=
                    settings.minimum_capture_fraction,
                minimum_utility=
                    settings.minimum_cost_aware_utility,
                condition_penalty_weight=
                    settings.condition_penalty_weight,
                pair_seconds_per_factorization=
                    rates["pair_seconds_per_factorization"],
                cayley_seconds_per_cubic_unit=
                    rates["cayley_seconds_per_cubic_unit"],
            )

        if not cost_ranked:
            events.append({
                "kind":"enrichment_rejected_v15",
                "step":int(step),
                "time":float(step*dt),
                "reason":"cost_aware_utility_gate",
                "basis_size":int(len(basis)),
                "relative_defect":
                    float(defect.relative_to_hpsi),
                "candidate_count":int(len(dynamic)),
                "residual_best_capture_fraction":
                    float(ranked[0].capture_fraction),
            })
            return False,None

        best_cost=cost_ranked[0]
        by_index={
            s.candidate_index:s
            for s in ranked
        }
        best=by_index[
            best_cost.candidate_index
        ]
        item=dynamic[
            best.candidate_index
        ]

        cost_history.append({
            "step":int(step),
            "time":float(step*dt),
            "candidate_index":
                int(best.candidate_index),
            "label":str(best.label),
            "capture_fraction":
                float(best.capture_fraction),
            "utility":
                float(best_cost.utility),
            "normalized_incremental_cost":
                float(
                    best_cost.normalized_incremental_cost
                ),
            "estimated_incremental_seconds":
                best_cost.estimated_incremental_seconds,
            "expanded_condition_number":
                float(best.expanded_condition_number),
            "residual_only_best_label":
                str(ranked[0].label),
            "residual_only_best_capture_fraction":
                float(ranked[0].capture_fraction),
        })

        child=item.candidate.to_tbf(
            next_uid,
            node_prefix="v15_cost_aware",
        )
        basis_before=len(basis)

        # Reuse the candidate's already-factorized expanded cache.
        expanded_cache=best.expanded_cache
        expanded_cache.basis[-1]=child

        # v0.14 would rebuild every old-old pair after insertion.
        ledger.v14_factorization_baseline += (
            v14_factorization_equivalent_for_sh(
                basis_before+1
            )
        )

        before=_cache_stat_snapshot(
            expanded_cache
        )
        with ledger.timed("matrix_build"):
            Snew,Hnew,Snucnew=(
                expand_cached_spinor_lvc_matrices(
                    S,H,Snuc,
                    expanded_cache,
                    provider,
                )
            )
        # All child pairs were already computed during candidate conditioning.
        ledger.add_cache_delta(
            expanded_cache.stats,
            before,
            category="propagation",
        )

        basis.append(child)
        C=np.concatenate([
            C,
            np.zeros(2,dtype=complex),
        ])
        S,H,Snuc,cache=(
            Snew,Hnew,Snucnew,expanded_cache
        )
        ledger.incremental_expansions+=1

        birth_step[int(next_uid)]=int(step)
        lineage[int(next_uid)]={
            "parent_uid":int(item.parent_uid),
            "birth_step":int(step),
            "birth_time":float(step*dt),
            "guidance_state":int(child.state),
            "source":"cost_aware_tdse_defect",
            "candidate_label":
                str(item.candidate.label),
            "predicted_capture_fraction":
                float(best.capture_fraction),
            "cost_aware_utility":
                float(best_cost.utility),
            "normalized_incremental_cost":
                float(
                    best_cost.normalized_incremental_cost
                ),
        }
        next_uid+=1

        renormalize()
        ledger.enrichment_events+=1
        last_adaptation_step=int(step)

        after=evaluate_defect(step)

        events.append({
            "kind":"cost_aware_defect_enrichment",
            "step":int(step),
            "time":float(step*dt),
            "parent_uid":int(item.parent_uid),
            "new_uid":int(child.uid),
            "guidance_state":int(child.state),
            "candidate_label":
                str(item.candidate.label),
            "basis_before":int(basis_before),
            "basis_after":int(len(basis)),
            "relative_defect_before":
                float(defect.relative_to_hpsi),
            "relative_defect_after":
                float(after.relative_to_hpsi),
            "defect_norm_before":
                float(defect.residual_norm),
            "defect_norm_after":
                float(after.residual_norm),
            "capture_fraction_predicted":
                float(best.capture_fraction),
            "residual_only_best_label":
                str(ranked[0].label),
            "residual_only_best_capture_fraction":
                float(ranked[0].capture_fraction),
            "cost_aware_utility":
                float(best_cost.utility),
            "normalized_incremental_cost":
                float(
                    best_cost.normalized_incremental_cost
                ),
            "estimated_incremental_seconds":
                best_cost.estimated_incremental_seconds,
            "expanded_condition_number":
                float(best.expanded_condition_number),
            "candidate_count":int(len(dynamic)),
            "residual_shortlist":
                int(len(ranked)),
            "zero_coefficient_insertion":True,
            "incremental_matrix_expansion":True,
            "new_pair_factorizations_during_expansion":
                int(
                    expanded_cache.stats.canonical_solves
                    -before["canonical_solves"]
                ),
        })
        return True,after

    def control(step):
        nonlocal latest_defect
        nonlocal low_defect_streak
        nonlocal last_adaptation_step

        latest_defect=evaluate_defect(step)

        if (
            current_condition()
            >settings.hard_condition_limit
            and len(basis)>settings.min_basis
        ):
            changed=pruning_attempt(
                step,
                "hard_condition_limit",
                settings.emergency_prune_fractional_loss,
                require_condition_improvement=True,
            )
            if changed:
                latest_defect=evaluate_defect(step)

        separated=(
            step-last_adaptation_step
            >=settings.minimum_adaptation_separation_steps
        )

        if (
            latest_defect.relative_to_hpsi
            >=settings.enrich_relative_threshold
        ):
            low_defect_streak=0
            if not separated:
                return

            if len(basis)>=settings.max_basis:
                changed=pruning_attempt(
                    step,
                    "basis_budget_replacement",
                    settings.max_replacement_prune_fractional_loss,
                    require_condition_improvement=False,
                )
                if changed:
                    latest_defect=evaluate_defect(step)

            if len(basis)<settings.max_basis:
                changed,after=enrichment_attempt(
                    step,latest_defect
                )
                if changed:
                    latest_defect=after

        elif (
            latest_defect.relative_to_hpsi
            <=settings.prune_relative_threshold
        ):
            low_defect_streak+=1
            if (
                low_defect_streak
                >=settings.prune_patience_checks
                and separated
                and len(basis)>settings.min_basis
            ):
                changed=pruning_attempt(
                    step,
                    "sustained_low_defect",
                    settings.max_prune_fractional_loss,
                    require_condition_improvement=False,
                )
                if changed:
                    low_defect_streak=0
                    latest_defect=evaluate_defect(step)
        else:
            low_defect_streak=0

    def record(step):
        rho=reduced_density_from_snuc(
            C,Snuc,normalize=True
        )
        records.append({
            "step":int(step),
            "time":float(step*dt),
            "norm":float(
                spinor_complete_generalized_norm(
                    C,S
                )
            ),
            "basis_size":int(len(basis)),
            "electronic_dimension":int(len(C)),
            "condition_number_nuclear":
                current_condition(),
            "diabatic_populations":
                np.real(np.diag(rho)),
            "diabatic_coherence":
                complex(rho[0,1]),
            "latest_relative_defect":
                None if latest_defect is None
                else float(
                    latest_defect.relative_to_hpsi
                ),
            "latest_defect_norm":
                None if latest_defect is None
                else float(
                    latest_defect.residual_norm
                ),
            "connection_correction_norm":
                float(latest_connection_correction),
            "connection_seed_norm":
                float(latest_connection_seed),
        })

    ledger.observe_basis(len(basis))

    if settings.check_initial_defect:
        control(0)
        ledger.observe_basis(len(basis))

    record(0)

    for step in range(1,int(steps)+1):
        old_basis=_copy_basis(basis)
        S0=S.copy()
        H0=H.copy()

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

        # New endpoint snapshot: one pair factorization per canonical pair, reused by
        # overlap, kinetic, and potential matrix elements.
        S1,H1,Snuc1,cache1=full_endpoint_build()

        mid=_midpoint_basis(
            old_basis,basis
        )
        qdot_mid=0.5*(qdot0+qdot1)
        pdot_mid=0.5*(pdot0+pdot1)

        mid_cache=GaussianPairCache(mid)
        ledger.pair_snapshots+=1
        ledger.v14_factorization_baseline += (
            v14_factorization_equivalent_for_time(
                len(mid)
            )
        )
        before=_cache_stat_snapshot(mid_cache)
        ledger.time_matrix_calls+=1
        with ledger.timed("time_matrix"):
            T_seed=build_cached_spinor_time_matrix(
                mid_cache,
                qdot_mid,
                pdot_mid,
            )
        ledger.add_cache_delta(
            mid_cache.stats,before,
            category="propagation",
        )

        T_mid=metric_compatible_basis_connection(
            S0,S1,dt,seed=T_seed
        )
        latest_connection_correction=float(
            np.linalg.norm(
                T_mid-T_seed,ord="fro"
            )
        )
        latest_connection_seed=float(
            np.linalg.norm(
                T_seed,ord="fro"
            )
        )

        ledger.cayley_solve_calls+=1
        ledger.cayley_cubic_units+=(
            dense_cubic_units(len(basis))
        )
        with ledger.timed("cayley_solve"):
            C=moving_basis_midpoint_cayley_step(
                C,S0,H0,S1,H1,T_mid,dt
            )

        S,H,Snuc,cache=(
            S1,H1,Snuc1,cache1
        )

        if step%settings.defect_interval==0:
            control(step)

        ledger.observe_basis(len(basis))
        basis_size_sum+=len(basis)
        basis_size_samples+=1

        if step%int(store_every)==0:
            record(step)

    if records[-1]["step"]!=int(steps):
        record(int(steps))

    ledger.stop()
    ledger.finalize_avoided()

    average_basis=float(
        basis_size_sum/max(
            basis_size_samples,1
        )
    )

    return {
        "records":records,
        "defect_history":defect_history,
        "cost_history":cost_history,
        "events":events,
        "final_basis":basis,
        "final_coefficients":C,
        "final_overlap":S.copy(),
        "final_hamiltonian":H.copy(),
        "final_nuclear_overlap":Snuc.copy(),
        "lineage":lineage,
        "average_basis_size":
            average_basis,
        "complexity":ledger.as_dict(),
        "settings":{
            "dt":float(dt),
            "steps":int(steps),
            "store_every":int(store_every),
            "control":asdict(settings),
            "representation":
                "time-adaptive spinor-complete global-diabatic Gaussian basis",
            "growth":
                "cost-aware TDSE-defect capture with energy-conserving local candidates",
            "matrix_algebra":
                "shared Gaussian-pair cache + incremental add/remove matrix updates",
            "pruning":
                "exact leave-one-out represented-wavefunction projection loss",
        },
    }
