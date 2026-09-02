import numpy as np

from .adaptive_spawning import CouplingExposureTracker
from .dynamic_gauge_graph import AnalyticCI2DFrameProvider
from .dynamic_graph_aims import DynamicGraphTBF, _verlet, _kinematics
from .managed_graph_aims_v11 import _spawn_from_integrated_exposure
from .moving_graph_gaussian import (
    metric_compatible_basis_connection,
    moving_basis_coefficient_step,
)
from .moving_basis_v12 import moving_basis_midpoint_cayley_step
from .spinor_complete_lvc_v12 import (
    build_spinor_complete_lvc_matrices,
    build_spinor_complete_time_matrix,
    spinor_complete_reduced_density,
    spinor_complete_generalized_norm,
    coefficients_matrix,
    flatten_coefficients,
)
from .paired_basis_management_v12 import (
    prune_nuclear_gaussian_pairs,
)


def _copy_basis(basis):
    return [
        DynamicGraphTBF(
            uid=int(b.uid),
            state=int(b.state),
            q=np.asarray(b.q,float).copy(),
            p=np.asarray(b.p,float).copy(),
            A=np.asarray(b.A,float).copy(),
            node=b.node,
            spawned_targets=set(b.spawned_targets),
        )
        for b in basis
    ]


def _midpoint_basis(old_basis,new_basis):
    if len(old_basis)!=len(new_basis):
        raise ValueError("midpoint basis requires equal sizes.")

    out=[]
    for old,new in zip(old_basis,new_basis):
        out.append(
            DynamicGraphTBF(
                uid=old.uid,
                state=old.state,
                q=0.5*(old.q+new.q),
                p=0.5*(old.p+new.p),
                A=old.A.copy(),
                node=old.node,
                spawned_targets=set(old.spawned_targets),
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
    return np.asarray(qdots),np.asarray(pdots)


def initialize_spinor_complete_coefficients(
    basis,
    provider,
    nuclear_coefficients=None,
):
    """Initialize each Gaussian with its local adiabatic guidance-state spinor."""
    n=len(basis)
    if nuclear_coefficients is None:
        nuclear_coefficients=np.ones(n,dtype=complex)
    nuclear_coefficients=np.asarray(nuclear_coefficients,dtype=complex)
    if nuclear_coefficients.shape!=(n,):
        raise ValueError("nuclear_coefficients must have one value per Gaussian.")

    C=np.zeros((n,2),dtype=complex)
    for i,b in enumerate(basis):
        point=provider.evaluate(np.asarray(b.q,float))
        C[i]=nuclear_coefficients[i]*point.frame[:,int(b.state)]
    return flatten_coefficients(C)


def run_spinor_complete_lvc_gaussians(
    initial_basis,
    C0=None,
    provider=None,
    dt=0.005,
    steps=120,
    integrator="cayley",
    spawn_action_threshold=1e-4,
    spawn_coupling_floor=1e-8,
    overlap_block=0.9999,
    child_overlap_block=0.995,
    max_basis=10,
    max_generation=5,
    children_per_event=2,
    allow_repeated_spawning=True,
    minimum_spawn_separation_steps=4,
    position_shifts=(0.0,0.05,-0.05),
    width_scales=(0.65,1.0,1.55),
    momentum_directions=("nac","momentum"),
    novelty_power=0.5,
    condition_limit=1e9,
    eigenvalue_floor=1e-10,
    max_pruning_loss=1e-7,
    store_every=20,
):
    r"""Spinor-complete Gaussian dynamics in the exact global diabatic LVC basis.

    The electronic subspace attached to every nuclear Gaussian is complete:

        Psi(R,t) = sum_k sum_a C_ka(t) g_k(R,t) |a_d>.

    Thus adiabatic/diabatic electronic rotations do not depend on whether a partner
    state happens to have been spawned at exactly the same nuclear phase-space point.

    The nuclear Gaussian center still carries a `state` used only for classical
    guidance and spawning logic.  The quantum electronic amplitudes are the full
    two-component coefficient block C_k.

    This is a classically guided vector-Gaussian benchmark model, not standard AIMS.
    """
    if integrator not in {"cayley","rk4"}:
        raise ValueError("integrator must be 'cayley' or 'rk4'.")

    provider=provider or AnalyticCI2DFrameProvider()
    if not isinstance(provider,AnalyticCI2DFrameProvider):
        raise TypeError(
            "The v0.12 spinor-complete exact-LVC runner is specialized to "
            "AnalyticCI2DFrameProvider."
        )

    basis=_copy_basis(initial_basis)
    if len(basis)==0:
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
            raise ValueError(
                "C0 must have shape (n,), (n,2), or (2*n,)."
            )

    S,H,Snuc=build_spinor_complete_lvc_matrices(
        basis,provider
    )
    C/=np.sqrt(spinor_complete_generalized_norm(C,S))

    lineage={
        int(b.uid):{
            "parent_uid":None,
            "generation":0,
            "birth_step":0,
            "birth_time":0.0,
            "birth_state":int(b.state),
            "width_scale_from_parent":1.0,
        }
        for b in basis
    }

    tracker=CouplingExposureTracker(
        action_threshold=spawn_action_threshold,
        coupling_floor=spawn_coupling_floor,
    )
    next_uid=max(int(b.uid) for b in basis)+1
    last_spawn_step={}
    events=[]
    records=[]

    candidate_kwargs=dict(
        position_shifts=position_shifts,
        width_scales=width_scales,
        momentum_directions=momentum_directions,
        overlap_block=overlap_block,
        child_overlap_block=child_overlap_block,
        novelty_power=novelty_power,
    )

    latest_connection_correction=0.0
    latest_connection_seed=0.0

    def record(step):
        rho=spinor_complete_reduced_density(
            C,basis,normalize=True
        )
        generations=[
            lineage[int(b.uid)]["generation"]
            for b in basis
        ]
        records.append({
            "step":int(step),
            "time":float(step*dt),
            "norm":float(spinor_complete_generalized_norm(C,S)),
            "basis_size":int(len(basis)),
            "electronic_dimension":int(len(C)),
            "condition_number_nuclear":float(np.linalg.cond(Snuc)),
            "maximum_generation":int(max(generations,default=0)),
            "diabatic_populations":np.real(np.diag(rho)),
            "diabatic_coherence":complex(rho[0,1]),
            "connection_correction_norm":float(latest_connection_correction),
            "connection_seed_norm":float(latest_connection_seed),
        })

    record(0)

    for step in range(1,int(steps)+1):
        old_basis=_copy_basis(basis)
        S0=S.copy()
        H0=H.copy()

        qdot0,pdot0=_kinematic_arrays(old_basis,provider)

        for b in basis:
            b.q,b.p=_verlet(b,provider,dt)

        qdot1,pdot1=_kinematic_arrays(basis,provider)

        S1,H1,Snuc1=build_spinor_complete_lvc_matrices(
            basis,provider
        )

        mid=_midpoint_basis(old_basis,basis)
        qdot_mid=0.5*(qdot0+qdot1)
        pdot_mid=0.5*(pdot0+pdot1)

        T_seed=build_spinor_complete_time_matrix(
            mid,qdot_mid,pdot_mid
        )
        T_mid=metric_compatible_basis_connection(
            S0,S1,dt,seed=T_seed
        )

        latest_connection_correction=np.linalg.norm(
            T_mid-T_seed,ord="fro"
        )
        latest_connection_seed=np.linalg.norm(
            T_seed,ord="fro"
        )

        if integrator=="cayley":
            C=moving_basis_midpoint_cayley_step(
                C,S0,H0,S1,H1,T_mid,dt
            )
        else:
            C=moving_basis_coefficient_step(
                C,S0,H0,S1,H1,T_mid,dt
            )

        S,H,Snuc=S1,H1,Snuc1

        if len(basis)<int(max_basis):
            spawn=_spawn_from_integrated_exposure(
                basis,
                provider,
                lineage,
                next_uid,
                tracker,
                dt,
                step,
                max_generation=max_generation,
                children_per_event=min(
                    int(children_per_event),
                    int(max_basis-len(basis)),
                ),
                last_spawn_step=last_spawn_step,
                minimum_spawn_separation_steps=minimum_spawn_separation_steps,
                allow_repeated_spawning=allow_repeated_spawning,
                candidate_kwargs=candidate_kwargs,
            )

            if spawn is not None:
                Cmat=coefficients_matrix(C,len(basis))

                for child,detail in zip(
                    spawn["children"],spawn["details"]
                ):
                    basis.append(child)
                    Cmat=np.vstack([
                        Cmat,
                        np.zeros((1,2),dtype=complex),
                    ])
                    lineage[int(child.uid)]={
                        "parent_uid":int(spawn["parent_uid"]),
                        "generation":int(detail["generation"]),
                        "birth_step":int(step),
                        "birth_time":float(step*dt),
                        "birth_state":int(child.state),
                        "width_scale_from_parent":float(detail["width_scale"]),
                    }

                C=flatten_coefficients(Cmat)
                next_uid=int(spawn["next_uid"])

                events.append({
                    "kind":"optimized_spawn",
                    "step":int(step),
                    "time":float(step*dt),
                    "parent_uid":int(spawn["parent_uid"]),
                    "parent_generation":int(spawn["parent_generation"]),
                    "target_state":int(spawn["target_state"]),
                    "integrated_coupling_action":float(
                        spawn["integrated_coupling_action"]
                    ),
                    "instantaneous_coupling_rate":float(
                        spawn["instantaneous_coupling_rate"]
                    ),
                    "children":spawn["details"],
                })

                S,H,Snuc=build_spinor_complete_lvc_matrices(
                    basis,provider
                )

        if len(basis)>1 and np.linalg.cond(Snuc)>condition_limit:
            Cmat=coefficients_matrix(C,len(basis))
            pruning=prune_nuclear_gaussian_pairs(
                Cmat,
                Snuc,
                condition_limit=condition_limit,
                eigenvalue_floor=eigenvalue_floor,
                max_projection_loss=max_pruning_loss,
            )

            if pruning.removed:
                removed=[basis[i].uid for i in pruning.removed]
                basis=[basis[i] for i in pruning.keep]
                C=flatten_coefficients(
                    pruning.coefficients_matrix
                )
                S,H,Snuc=build_spinor_complete_lvc_matrices(
                    basis,provider
                )

                events.append({
                    "kind":"paired_prune",
                    "step":int(step),
                    "time":float(step*dt),
                    "removed_uids":[int(x) for x in removed],
                    "projection_loss":float(pruning.projection_loss),
                    "condition_before":float(pruning.condition_before),
                    "condition_after":float(pruning.condition_after),
                })

        if step%int(store_every)==0:
            record(step)

    if records[-1]["step"]!=int(steps):
        record(int(steps))

    return {
        "records":records,
        "events":events,
        "final_basis":basis,
        "final_coefficients":C,
        "final_overlap":S.copy(),
        "final_hamiltonian":H.copy(),
        "final_nuclear_overlap":Snuc.copy(),
        "lineage":lineage,
        "settings":{
            "dt":float(dt),
            "steps":int(steps),
            "integrator":integrator,
            "representation":"spinor-complete global diabatic electronic basis",
            "max_basis":int(max_basis),
            "max_generation":int(max_generation),
            "children_per_event":int(children_per_event),
            "spawn_action_threshold":float(spawn_action_threshold),
            "position_shifts":[float(x) for x in position_shifts],
            "width_scales":[float(x) for x in width_scales],
        },
    }
