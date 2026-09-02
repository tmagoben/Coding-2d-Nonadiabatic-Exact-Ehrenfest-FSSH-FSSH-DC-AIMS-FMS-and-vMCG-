from dataclasses import dataclass, asdict
import numpy as np

from .dynamic_gauge_graph import AnalyticCI2DFrameProvider
from .dynamic_graph_aims import DynamicGraphTBF, _verlet, _kinematics
from .moving_graph_gaussian import metric_compatible_basis_connection
from .moving_basis_v12 import moving_basis_midpoint_cayley_step
from .spinor_complete_lvc_v12 import (
    build_spinor_complete_time_matrix,
    spinor_complete_reduced_density,
    spinor_complete_generalized_norm,
    coefficients_matrix,
    flatten_coefficients,
)
from .spinor_complete_dynamics_v12 import initialize_spinor_complete_coefficients
from .tdse_defect_v13 import (
    TDSEDefect,
    reconstruct_spinor_complete_wavefunction,
    reconstruct_spinor_complete_time_derivative,
    apply_lvc_grid_hamiltonian,
)
from .gaussian_nd import gaussian_nd
from .residual_basis_v13 import nuclear_overlap_matrix
from .defect_candidates_v14 import (
    generate_energy_conserving_defect_candidates,
    rank_dynamic_defect_candidates_prepared,
)
from .residual_pruning_v14 import prune_low_loss_gaussian_pair
from .complexity_v14 import ComplexityLedger
from .fast_lvc_matrices_v14 import (
    build_spinor_complete_lvc_matrices_symmetric,
    hermitian_pair_evaluation_count,
)


@dataclass(frozen=True)
class AdaptiveDefectSettings:
    """Error-control policy for the v0.14 adaptive reference implementation."""

    defect_interval: int = 10
    enrich_relative_threshold: float = 0.025
    prune_relative_threshold: float = 0.010
    minimum_capture_fraction: float = 0.015

    min_basis: int = 4
    max_basis: int = 12

    minimum_adaptation_separation_steps: int = 10
    minimum_prune_age_steps: int = 20
    prune_patience_checks: int = 2

    max_prune_fractional_loss: float = 2e-6
    max_replacement_prune_fractional_loss: float = 2e-6
    emergency_prune_fractional_loss: float = 1e-4

    condition_limit: float = 2e5
    hard_condition_limit: float = 1e7
    orthogonal_norm_floor: float = 1e-8

    candidate_position_shifts: tuple = (0.0, 0.06, -0.06)
    candidate_width_scales: tuple = (0.75, 1.0, 1.35)
    candidate_momentum_directions: tuple = ("nac", "momentum")
    include_same_surface_candidates: bool = True
    include_other_surface_candidates: bool = True
    candidate_overlap_block: float = 0.999999

    check_initial_defect: bool = True

    def validate(self):
        if self.defect_interval<=0:
            raise ValueError("defect_interval must be positive.")
        if self.enrich_relative_threshold<=self.prune_relative_threshold:
            raise ValueError(
                "enrichment threshold must exceed pruning threshold to create "
                "hysteresis."
            )
        if self.min_basis<1 or self.max_basis<self.min_basis:
            raise ValueError("invalid min/max basis.")
        if not (0.0<=self.minimum_capture_fraction<=1.0):
            raise ValueError("minimum_capture_fraction must be in [0,1].")
        if self.condition_limit<=1.0 or self.hard_condition_limit<self.condition_limit:
            raise ValueError("invalid condition-number limits.")
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
            spawned_targets=set(getattr(b,"spawned_targets",set())),
        )
        for b in basis
    ]


def _midpoint_basis(old_basis,new_basis):
    out=[]
    if len(old_basis)!=len(new_basis):
        raise ValueError("midpoint basis requires equal basis sizes.")
    for old,new in zip(old_basis,new_basis):
        if old.uid!=new.uid:
            raise ValueError("basis identity changed inside one propagation step.")
        out.append(
            DynamicGraphTBF(
                uid=old.uid,
                state=old.state,
                q=0.5*(old.q+new.q),
                p=0.5*(old.p+new.p),
                A=old.A.copy(),
                node=old.node,
                spawned_targets=set(getattr(old,"spawned_targets",set())),
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


def compute_tdse_defect_with_matrices(
    C,
    basis,
    provider,
    grid,
    S,
    H,
    Snuc=None,
):
    r"""Evaluate the v0.13 TDSE defect while reusing already-built S and H.

    This avoids rebuilding the pairwise exact-LVC matrices at every defect checkpoint.
    The additional work is:

    - one moving-basis T matrix;
    - one dense projected linear solve;
    - wavefunction/Psidot reconstruction on the diagnostic grid;
    - one FFT kinetic application.
    """
    C=np.asarray(C,dtype=complex)
    if Snuc is None:
        Snuc=nuclear_overlap_matrix(basis)

    qdots,pdots=_kinematic_arrays(basis,provider)
    T=build_spinor_complete_time_matrix(
        basis,qdots,pdots
    )
    Cdot=np.linalg.solve(
        np.asarray(S,dtype=complex),
        -(1j*np.asarray(H,dtype=complex)+T)@C,
    )

    psi=reconstruct_spinor_complete_wavefunction(
        C,basis,grid.points
    )
    psidot=reconstruct_spinor_complete_time_derivative(
        C,Cdot,basis,grid.points,qdots,pdots
    )
    Hpsi=apply_lvc_grid_hamiltonian(psi,grid)
    residual=1j*psidot-Hpsi

    r2=float(max(
        np.real(_grid_inner(residual,residual,grid.area)),
        0.0,
    ))
    h2=float(max(
        np.real(_grid_inner(Hpsi,Hpsi,grid.area)),
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
        projected += float(np.real(
            np.vdot(b[:,a],coeff)
        ))

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


def run_time_adaptive_defect_lvc_gaussians(
    initial_basis,
    C0=None,
    provider=None,
    grid=None,
    dt=0.005,
    steps=120,
    settings=AdaptiveDefectSettings(),
    store_every=10,
):
    r"""Fully time-adaptive residual-controlled Gaussian propagation.

    Control loop
    ------------
    1. Propagate the current moving spinor-complete Gaussian basis with the
       generalized midpoint/Cayley step.
    2. Every `defect_interval` steps evaluate

           R = i dPsi/dt - H Psi

       on the exact diagnostic grid.
    3. If the relative defect exceeds the enrichment threshold, generate a local,
       energy-conserving candidate dictionary and select the Gaussian that captures
       the largest unresolved defect.
    4. Insert the new Gaussian with a zero two-component electronic coefficient.
    5. If the defect remains small for multiple checks, remove a Gaussian only when
       its exact leave-one-out wavefunction projection loss is below the pruning
       budget.
    6. Use separate enrichment/pruning thresholds plus an adaptation cooldown to
       prevent basis-size thrashing.

    Scientific scope
    ----------------
    This is a low-dimensional analytic-LVC reference algorithm.  It is not a
    production molecular AIMS implementation.
    """
    settings=settings.validate()
    provider=provider or AnalyticCI2DFrameProvider()

    if not isinstance(provider,AnalyticCI2DFrameProvider):
        raise TypeError(
            "v0.14 adaptive defect propagation is specialized to the analytic "
            "two-state LVC provider."
        )
    if grid is None:
        raise ValueError("A diagnostic grid is required for TDSE-defect control.")
    if int(store_every)<=0:
        raise ValueError("store_every must be positive.")

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

    ledger=ComplexityLedger()
    ledger.start()

    def build_matrices():
        ledger.matrix_build_calls+=1
        ledger.pair_matrix_evaluations += hermitian_pair_evaluation_count(
            len(basis)
        )
        ledger.ordered_pair_equivalent += int(len(basis)**2)
        with ledger.timed("matrix_build"):
            return build_spinor_complete_lvc_matrices_symmetric(
                basis,provider
            )

    S,H,Snuc=build_matrices()
    C/=np.sqrt(
        spinor_complete_generalized_norm(C,S)
    )

    next_uid=max(int(b.uid) for b in basis)+1
    birth_step={int(b.uid):0 for b in basis}
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

    latest_defect=None
    latest_defect_after=None
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
            if step-int(born)<settings.minimum_prune_age_steps
        }

    def renormalize():
        nonlocal C
        norm=spinor_complete_generalized_norm(C,S)
        if norm<=0.0 or not np.isfinite(norm):
            raise RuntimeError("invalid generalized norm.")
        C=C/np.sqrt(norm)

    def pruning_attempt(step, reason, max_loss, require_condition_improvement=False):
        nonlocal basis,C,S,H,Snuc,last_adaptation_step
        if len(basis)<=settings.min_basis:
            return False

        ledger.pruning_audits+=1
        with ledger.timed("pruning"):
            Cmat=coefficients_matrix(C,len(basis))
            result=prune_low_loss_gaussian_pair(
                Cmat,
                Snuc,
                uids=[b.uid for b in basis],
                max_fractional_loss=max_loss,
                protected_uids=protected_uids(step),
                require_condition_improvement=require_condition_improvement,
            )

        if result is None:
            return False

        old_basis_size=len(basis)
        removed_uid=int(result.removed_uid)
        basis=[basis[i] for i in result.keep]
        C=flatten_coefficients(
            result.coefficients_matrix
        )
        S,H,Snuc=build_matrices()
        renormalize()

        birth_step.pop(removed_uid,None)
        lineage.setdefault(removed_uid,{})["removed_step"]=int(step)
        lineage[removed_uid]["removed_time"]=float(step*dt)
        lineage[removed_uid]["removal_reason"]=str(reason)

        ledger.pruning_events+=1
        last_adaptation_step=int(step)

        events.append({
            "kind":"residual_prune",
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
        })
        return True

    def evaluate_defect(step):
        ledger.defect_evaluations+=1
        with ledger.timed("defect"):
            d=compute_tdse_defect_with_matrices(
                C,basis,provider,grid,S,H,Snuc
            )
        defect_history.append({
            "step":int(step),
            "time":float(step*dt),
            "basis_size":int(len(basis)),
            "residual_norm":float(d.residual_norm),
            "relative_to_hpsi":float(d.relative_to_hpsi),
            "projected_residual_norm":
                float(d.projected_residual_norm),
            "condition_number":current_condition(),
        })
        return d

    def enrichment_attempt(step, defect):
        nonlocal basis,C,S,H,Snuc,next_uid,last_adaptation_step

        dynamic=generate_energy_conserving_defect_candidates(
            basis,
            provider,
            position_shifts=settings.candidate_position_shifts,
            width_scales=settings.candidate_width_scales,
            momentum_directions=settings.candidate_momentum_directions,
            include_same_surface=
                settings.include_same_surface_candidates,
            include_other_surfaces=
                settings.include_other_surface_candidates,
            overlap_block=settings.candidate_overlap_block,
        )

        ledger.candidate_ranking_calls+=1
        ledger.candidate_count_scored+=len(dynamic)
        ledger.peak_candidate_count=max(
            ledger.peak_candidate_count,len(dynamic)
        )

        with ledger.timed("candidate_ranking"):
            ranked=rank_dynamic_defect_candidates_prepared(
                defect,
                basis,
                dynamic,
                grid,
                condition_limit=settings.condition_limit,
                orthogonal_norm_floor=
                    settings.orthogonal_norm_floor,
                max_return=4,
            )

        if not ranked:
            events.append({
                "kind":"enrichment_rejected",
                "step":int(step),
                "time":float(step*dt),
                "reason":"no_admissible_candidate",
                "basis_size":int(len(basis)),
                "relative_defect":float(defect.relative_to_hpsi),
                "candidate_count":int(len(dynamic)),
            })
            return False,None

        best=ranked[0]
        if best.capture_fraction<settings.minimum_capture_fraction:
            events.append({
                "kind":"enrichment_rejected",
                "step":int(step),
                "time":float(step*dt),
                "reason":"insufficient_capture_fraction",
                "basis_size":int(len(basis)),
                "relative_defect":float(defect.relative_to_hpsi),
                "capture_fraction":float(best.capture_fraction),
                "candidate_count":int(len(dynamic)),
            })
            return False,None

        item=dynamic[best.candidate_index]
        child=item.candidate.to_tbf(
            next_uid,
            node_prefix="v14_defect",
        )
        basis_before=len(basis)
        basis.append(child)
        C=np.concatenate([
            C,
            np.zeros(2,dtype=complex),
        ])

        birth_step[int(next_uid)]=int(step)
        lineage[int(next_uid)]={
            "parent_uid":int(item.parent_uid),
            "birth_step":int(step),
            "birth_time":float(step*dt),
            "guidance_state":int(child.state),
            "source":"tdse_defect",
            "candidate_label":str(item.candidate.label),
            "predicted_capture_fraction":
                float(best.capture_fraction),
        }
        next_uid+=1

        S,H,Snuc=build_matrices()
        renormalize()

        ledger.enrichment_events+=1
        last_adaptation_step=int(step)

        after=evaluate_defect(step)

        events.append({
            "kind":"defect_enrichment",
            "step":int(step),
            "time":float(step*dt),
            "parent_uid":int(item.parent_uid),
            "new_uid":int(child.uid),
            "guidance_state":int(child.state),
            "candidate_label":str(item.candidate.label),
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
            "captured_defect_norm_predicted":
                float(best.captured_defect_norm),
            "orthogonal_norm":
                float(best.orthogonal_norm),
            "expanded_condition_number":
                float(best.expanded_condition_number),
            "candidate_count":int(len(dynamic)),
            "zero_coefficient_insertion":True,
        })
        return True,after

    def control(step):
        nonlocal latest_defect,latest_defect_after
        nonlocal low_defect_streak,last_adaptation_step

        latest_defect=evaluate_defect(step)
        latest_defect_after=None

        # Emergency conditioning control is independent of the residual thresholds.
        if (
            current_condition()>settings.hard_condition_limit
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

        if latest_defect.relative_to_hpsi>=settings.enrich_relative_threshold:
            low_defect_streak=0
            if not separated:
                return

            # At the hard basis budget, attempt a low-loss replacement first.
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
                    latest_defect_after=after

        elif latest_defect.relative_to_hpsi<=settings.prune_relative_threshold:
            low_defect_streak+=1
            if (
                low_defect_streak>=settings.prune_patience_checks
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
        rho=spinor_complete_reduced_density(
            C,basis,normalize=True
        )
        records.append({
            "step":int(step),
            "time":float(step*dt),
            "norm":float(
                spinor_complete_generalized_norm(C,S)
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
                else float(latest_defect.relative_to_hpsi),
            "latest_defect_norm":
                None if latest_defect is None
                else float(latest_defect.residual_norm),
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

        S1,H1,Snuc1=build_matrices()

        mid=_midpoint_basis(
            old_basis,basis
        )
        qdot_mid=0.5*(qdot0+qdot1)
        pdot_mid=0.5*(pdot0+pdot1)

        ledger.time_matrix_calls+=1
        with ledger.timed("time_matrix"):
            T_seed=build_spinor_complete_time_matrix(
                mid,qdot_mid,pdot_mid
            )
        T_mid=metric_compatible_basis_connection(
            S0,S1,dt,seed=T_seed
        )
        latest_connection_correction=float(
            np.linalg.norm(T_mid-T_seed,ord="fro")
        )
        latest_connection_seed=float(
            np.linalg.norm(T_seed,ord="fro")
        )

        ledger.cayley_solve_calls+=1
        with ledger.timed("cayley_solve"):
            C=moving_basis_midpoint_cayley_step(
                C,S0,H0,S1,H1,T_mid,dt
            )

        S,H,Snuc=S1,H1,Snuc1

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

    average_basis=float(
        basis_size_sum/max(basis_size_samples,1)
    )

    return {
        "records":records,
        "defect_history":defect_history,
        "events":events,
        "final_basis":basis,
        "final_coefficients":C,
        "final_overlap":S.copy(),
        "final_hamiltonian":H.copy(),
        "final_nuclear_overlap":Snuc.copy(),
        "lineage":lineage,
        "average_basis_size":average_basis,
        "complexity":ledger.as_dict(),
        "settings":{
            "dt":float(dt),
            "steps":int(steps),
            "store_every":int(store_every),
            "control":asdict(settings),
            "representation":
                "time-adaptive spinor-complete global-diabatic Gaussian basis",
            "growth":
                "TDSE-defect residual capture with energy-conserving local candidates",
            "pruning":
                "exact leave-one-out represented-wavefunction projection loss",
        },
    }
