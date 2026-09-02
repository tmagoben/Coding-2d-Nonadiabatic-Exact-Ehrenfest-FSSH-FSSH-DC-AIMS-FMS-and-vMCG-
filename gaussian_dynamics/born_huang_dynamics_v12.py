import numpy as np

from .adaptive_spawning import CouplingExposureTracker
from .basis_management import prune_redundant_basis
from .dynamic_gauge_graph import AnalyticCI2DFrameProvider
from .dynamic_graph_aims import DynamicGraphTBF, _verlet, _kinematics
from .managed_graph_aims_v11 import _spawn_from_integrated_exposure
from .moving_graph_gaussian import (
    metric_compatible_basis_connection,
    moving_basis_coefficient_step,
    generalized_norm,
)
from .moving_basis_v12 import moving_basis_midpoint_cayley_step
from .born_huang_grid_v12 import (
    build_born_huang_grid_2d,
    build_born_huang_matrices,
    born_huang_basis_time_matrix_grid,
    born_huang_reduced_density,
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
        raise ValueError("midpoint basis requires equal basis sizes.")

    out=[]
    for old,new in zip(old_basis,new_basis):
        if old.uid!=new.uid or old.state!=new.state:
            raise ValueError("basis identity changed inside one step.")
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


def run_born_huang_projected_gaussians(
    initial_basis,
    C0,
    provider=None,
    grid=None,
    grid_n=64,
    half_width=4.0,
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
    r"""Projected Born-Huang Gaussian dynamics for the analytic 2D CI benchmark.

    Basis:
        Xi_i(R) = g_i(R) Phi_ai(R),

    with the *coordinate-dependent* adiabatic electronic state over the full nuclear
    coordinate R.

    The basis Hamiltonian is evaluated by applying the same global-diabatic FFT-grid
    Hamiltonian used by the exact reference to every Xi_i and projecting back onto the
    Gaussian basis.  Therefore no center-frozen electronic approximation is made in
    this benchmark path.

    This is the most faithful v0.12 bridge between the adaptive Gaussian basis and the
    exact 2D TDSE.  It is a benchmark/reference implementation; its grid-projected
    matrix evaluation is not intended as the scalable molecular AIMS algorithm.
    """
    if integrator not in {"cayley","rk4"}:
        raise ValueError("integrator must be 'cayley' or 'rk4'.")

    provider=provider or AnalyticCI2DFrameProvider()
    if not isinstance(provider,AnalyticCI2DFrameProvider):
        raise TypeError(
            "The Born-Huang grid runner is specialized to AnalyticCI2DFrameProvider."
        )

    if grid is None:
        grid=build_born_huang_grid_2d(
            grid_n=grid_n,
            half_width=half_width,
            mass=provider.nuclear_mass_au,
            params=provider.params,
        )

    basis=_copy_basis(initial_basis)
    C=np.asarray(C0,dtype=complex).copy()
    if C.shape!=(len(basis),) or len(basis)==0:
        raise ValueError("C0 must contain one scalar coefficient per TBF.")

    S,H=build_born_huang_matrices(basis,grid)
    C/=np.sqrt(generalized_norm(C,S))

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
        rho=born_huang_reduced_density(
            C,basis,grid,normalize=True
        )
        generations=[
            lineage[int(b.uid)]["generation"]
            for b in basis
        ]

        records.append({
            "step":int(step),
            "time":float(step*dt),
            "norm":float(generalized_norm(C,S)),
            "basis_size":int(len(basis)),
            "condition_number":float(np.linalg.cond(S)),
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

        S1,H1=build_born_huang_matrices(basis,grid)

        mid=_midpoint_basis(old_basis,basis)
        qdot_mid=0.5*(qdot0+qdot1)
        pdot_mid=0.5*(pdot0+pdot1)

        T_seed=born_huang_basis_time_matrix_grid(
            mid,
            grid,
            qdot_mid,
            pdot_mid,
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

        S,H=S1,H1

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
                for child,detail in zip(
                    spawn["children"],spawn["details"]
                ):
                    basis.append(child)
                    C=np.concatenate([C,[0.0+0.0j]])
                    lineage[int(child.uid)]={
                        "parent_uid":int(spawn["parent_uid"]),
                        "generation":int(detail["generation"]),
                        "birth_step":int(step),
                        "birth_time":float(step*dt),
                        "birth_state":int(child.state),
                        "width_scale_from_parent":float(detail["width_scale"]),
                    }

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

                S,H=build_born_huang_matrices(basis,grid)

        if len(basis)>1 and np.linalg.cond(S)>condition_limit:
            pruning=prune_redundant_basis(
                C,S,
                condition_limit=condition_limit,
                eigenvalue_floor=eigenvalue_floor,
                max_projection_loss=max_pruning_loss,
            )

            if pruning.removed:
                removed=[basis[i].uid for i in pruning.removed]
                basis=[basis[i] for i in pruning.keep]
                C=pruning.coefficients
                S,H=build_born_huang_matrices(basis,grid)

                events.append({
                    "kind":"prune",
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
        "lineage":lineage,
        "grid":grid,
        "settings":{
            "dt":float(dt),
            "steps":int(steps),
            "integrator":integrator,
            "representation":"coordinate-dependent Born-Huang adiabatic TBFs",
            "hamiltonian_projection":"global-diabatic FFT-grid projection",
            "grid_n":int(grid.grid_n),
            "half_width":float(0.5*grid.grid_n*grid.dx),
            "max_basis":int(max_basis),
            "max_generation":int(max_generation),
            "children_per_event":int(children_per_event),
            "spawn_action_threshold":float(spawn_action_threshold),
        },
    }
