import numpy as np

from .adaptive_spawning import CouplingExposureTracker
from .basis_management import prune_redundant_basis
from .dynamic_gauge_graph import IncrementalElectronicGraph, AnalyticCI2DFrameProvider
from .dynamic_graph_aims import (
    DynamicGraphTBF,
    _center_node,
    _centroid_node,
    _verlet,
    _kinematics,
    _indicator,
)
from .graph_gaussian import GraphGaussianTBF
from .moving_graph_gaussian import (
    metric_compatible_basis_connection,
    moving_basis_coefficient_step,
    generalized_norm,
)
from .moving_graph_gaussian_v11 import nuclear_seed_basis_time_matrix_general
from .spa_matrix_elements_v11 import (
    build_graph_gaussian_matrices_spa_general,
    spa1_correction_norm_general,
)
from .optimized_spawning import select_spawn_children


def _ensure_graph(step,basis,manager,provider):
    for b in basis:
        node=_center_node(step,b.uid)
        if node not in manager.frames:
            connect=[b.node] if b.node in manager.frames and b.node!=node else []
            manager.add_from_provider(node,b.q,provider,connect_to=connect)
        b.node=node

    refs={}
    for i in range(len(basis)):
        refs[(i,i)]=basis[i].node
        for j in range(i+1,len(basis)):
            node=_centroid_node(step,basis[i].uid,basis[j].uid)
            if node not in manager.frames:
                manager.add_from_provider(
                    node,
                    0.5*(basis[i].q+basis[j].q),
                    provider,
                    connect_to=[basis[i].node,basis[j].node],
                )
            refs[(i,j)]=node
            refs[(j,i)]=node
    return refs


def _graph_basis(basis,dimension):
    out=[]
    for b in basis:
        c=np.zeros(dimension,dtype=complex)
        c[b.state]=1.0
        out.append(
            GraphGaussianTBF(
                b.node,
                b.q.copy(),
                b.p.copy(),
                b.A.copy(),
                c,
            )
        )
    return out


def _copy_basis(basis):
    return [
        DynamicGraphTBF(
            b.uid,
            b.state,
            b.q.copy(),
            b.p.copy(),
            b.A.copy(),
            b.node,
            set(b.spawned_targets),
        )
        for b in basis
    ]


def _spawn_from_integrated_exposure(
    basis,
    provider,
    lineage,
    next_uid,
    tracker,
    dt,
    current_step,
    max_generation,
    children_per_event,
    last_spawn_step,
    minimum_spawn_separation_steps,
    allow_repeated_spawning,
    candidate_kwargs,
):
    ns=provider.evaluate(basis[0].q).frame.shape[1]

    for parent in basis:
        generation=int(lineage[parent.uid]["generation"])
        if generation>=int(max_generation):
            continue

        for target in range(ns):
            if target==parent.state:
                continue
            if (not allow_repeated_spawning) and target in parent.spawned_targets:
                continue

            if allow_repeated_spawning:
                previous=last_spawn_step.get((parent.uid,target))
                if previous is not None:
                    if current_step-previous<int(minimum_spawn_separation_steps):
                        continue

            rate=_indicator(parent,target,provider)
            key=(int(parent.uid),int(target))
            ready,action=tracker.update(key,rate,dt)

            if not ready:
                continue

            selected=select_spawn_children(
                parent,
                target,
                provider,
                basis,
                children_per_event=children_per_event,
                **candidate_kwargs,
            )
            if not selected:
                continue

            tracker.consume(key)
            last_spawn_step[key]=int(current_step)

            if not allow_repeated_spawning:
                parent.spawned_targets.add(target)

            children=[]
            details=[]
            uid=next_uid

            for candidate in selected:
                child=DynamicGraphTBF(
                    uid=uid,
                    state=target,
                    q=candidate.q.copy(),
                    p=candidate.p.copy(),
                    A=candidate.A.copy(),
                    node=parent.node,
                )
                children.append(child)
                details.append({
                    "child_uid":int(uid),
                    "target_state":int(target),
                    "generation":generation+1,
                    "position_direction":candidate.position_direction,
                    "position_shift":candidate.position_shift,
                    "momentum_direction":candidate.momentum_direction,
                    "width_scale":candidate.width_scale,
                    "coupling_proxy":candidate.coupling_proxy,
                    "nuclear_overlap":candidate.nuclear_overlap,
                    "max_existing_overlap":candidate.max_existing_overlap,
                    "novelty":candidate.novelty,
                    "spawn_score":candidate.score,
                    "energy_residual":candidate.energy_residual,
                })
                uid+=1

            return {
                "parent_uid":int(parent.uid),
                "parent_generation":generation,
                "target_state":int(target),
                "integrated_coupling_action":float(action),
                "instantaneous_coupling_rate":float(rate),
                "children":children,
                "details":details,
                "next_uid":uid,
            }

    return None


def run_basis_complete_graph_aims(
    initial_basis,
    C0,
    provider=None,
    dt=0.0025,
    steps=100,
    spa_order=1,
    spawn_action_threshold=1e-4,
    spawn_coupling_floor=1e-8,
    overlap_block=0.9995,
    child_overlap_block=0.985,
    max_basis=12,
    max_generation=4,
    children_per_event=2,
    allow_repeated_spawning=True,
    minimum_spawn_separation_steps=4,
    position_shifts=(0.0,0.04,-0.04,0.08,-0.08),
    width_scales=(0.65,1.0,1.55),
    momentum_directions=("nac","momentum"),
    novelty_power=0.5,
    condition_limit=1e9,
    eigenvalue_floor=1e-10,
    max_pruning_loss=1e-7,
    store_every=10,
):
    """v0.11 basis-completeness graph-AIMS-style prototype.

    Differences from v0.10
    ----------------------
    - unequal-width Gaussian matrix elements are exact within SPA0/SPA1;
    - spawned children are selected from a local constrained phase-space search;
    - one event may add multiple nonredundant child candidates;
    - descendants may themselves spawn up to `max_generation`;
    - width diversity is a controlled basis-enrichment coordinate;
    - complete lineage and generation diagnostics are returned.

    This remains an analytic benchmark/prototype, not a production implementation of
    the Yang-Coe-Kaduk-Martinez optimal-spawning algorithm.
    """
    if spa_order not in (0,1):
        raise ValueError("spa_order must be 0 or 1.")
    if children_per_event<=0:
        raise ValueError("children_per_event must be positive.")
    if max_generation<0:
        raise ValueError("max_generation must be nonnegative.")

    provider=provider or AnalyticCI2DFrameProvider()
    dimension=provider.evaluate(initial_basis[0].q).frame.shape[1]
    manager=IncrementalElectronicGraph(dimension)

    basis=_copy_basis(initial_basis)
    C=np.asarray(C0,dtype=complex).copy()
    if len(C)!=len(basis):
        raise ValueError("C0 length mismatch.")

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

    refs=_ensure_graph(0,basis,manager,provider)
    gbasis=_graph_basis(basis,dimension)
    mass_matrix=provider.evaluate(basis[0].q).mass_matrix

    def build(order):
        return build_graph_gaussian_matrices_spa_general(
            _graph_basis(basis,dimension),
            manager.registry,
            mass_matrix,
            reference_selector=lambda i,j:refs[(i,j)],
            order=order,
        )

    S,H=build(spa_order)

    events=[]

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
            refs=_ensure_graph(0,basis,manager,provider)
            S,H=build(spa_order)
            events.append({
                "kind":"prune",
                "step":0,
                "time":0.0,
                "removed_uids":[int(x) for x in removed],
                "projection_loss":float(pruning.projection_loss),
                "condition_before":float(pruning.condition_before),
                "condition_after":float(pruning.condition_after),
            })

    C/=np.sqrt(generalized_norm(C,S))

    tracker=CouplingExposureTracker(
        action_threshold=spawn_action_threshold,
        coupling_floor=spawn_coupling_floor,
    )
    next_uid=max(b.uid for b in basis)+1
    last_spawn_step={}
    records=[]

    candidate_kwargs=dict(
        position_shifts=position_shifts,
        width_scales=width_scales,
        momentum_directions=momentum_directions,
        overlap_block=overlap_block,
        child_overlap_block=child_overlap_block,
        novelty_power=novelty_power,
    )

    def record(step):
        norm=generalized_norm(C,S)

        state_proxy=np.zeros(dimension)
        for state in range(dimension):
            idx=[i for i,b in enumerate(basis) if b.state==state]
            if idx:
                block=S[np.ix_(idx,idx)]
                cc=C[idx]
                state_proxy[state]=np.real(np.vdot(cc,block@cc))
        if norm>0:
            state_proxy/=norm

        H0=build_graph_gaussian_matrices_spa_general(
            _graph_basis(basis,dimension),
            manager.registry,
            mass_matrix,
            reference_selector=lambda i,j:refs[(i,j)],
            order=0,
        )[1]
        H1=build_graph_gaussian_matrices_spa_general(
            _graph_basis(basis,dimension),
            manager.registry,
            mass_matrix,
            reference_selector=lambda i,j:refs[(i,j)],
            order=1,
        )[1]

        generations=[lineage[b.uid]["generation"] for b in basis]
        width_dets=[float(np.linalg.det(b.A)) for b in basis]

        records.append({
            "step":int(step),
            "time":float(step*dt),
            "norm":float(norm),
            "basis_size":int(len(basis)),
            "state_population_proxy":state_proxy,
            "condition_number":float(np.linalg.cond(S)),
            "spa_order":int(spa_order),
            "spa1_relative_correction":spa1_correction_norm_general(H0,H1),
            "maximum_generation":int(max(generations,default=0)),
            "width_determinants":width_dets,
            "graph":manager.summary(),
        })

    record(0)

    for step in range(1,int(steps)+1):
        old=_copy_basis(basis)
        old_refs=refs
        old_gbasis=_graph_basis(old,dimension)
        S0,H0=S,H

        for b in basis:
            b.q,b.p=_verlet(b,provider,dt)

        refs=_ensure_graph(step,basis,manager,provider)
        S1,H1=build(spa_order)

        qdots=[]
        pdots=[]
        for b in old:
            qdot,pdot=_kinematics(b,provider)
            qdots.append(qdot)
            pdots.append(pdot)

        seed=nuclear_seed_basis_time_matrix_general(
            old_gbasis,
            manager.registry,
            lambda i,j:old_refs[(i,j)],
            np.asarray(qdots),
            np.asarray(pdots),
        )

        T=metric_compatible_basis_connection(
            S0,S1,dt,seed=seed
        )
        C=moving_basis_coefficient_step(
            C,S0,H0,S1,H1,T,dt
        )
        S,H=S1,H1

        if len(basis)<max_basis:
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
                for child,detail in zip(spawn["children"],spawn["details"]):
                    basis.append(child)
                    C=np.concatenate([C,[0.0+0.0j]])
                    lineage[child.uid]={
                        "parent_uid":spawn["parent_uid"],
                        "generation":detail["generation"],
                        "birth_step":int(step),
                        "birth_time":float(step*dt),
                        "birth_state":int(child.state),
                        "width_scale_from_parent":float(detail["width_scale"]),
                    }

                next_uid=spawn["next_uid"]

                events.append({
                    "kind":"optimized_spawn",
                    "step":int(step),
                    "time":float(step*dt),
                    "parent_uid":spawn["parent_uid"],
                    "parent_generation":spawn["parent_generation"],
                    "target_state":spawn["target_state"],
                    "integrated_coupling_action":spawn["integrated_coupling_action"],
                    "instantaneous_coupling_rate":spawn["instantaneous_coupling_rate"],
                    "children":spawn["details"],
                })

                refs=_ensure_graph(step,basis,manager,provider)
                S,H=build(spa_order)

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
                refs=_ensure_graph(step,basis,manager,provider)
                S,H=build(spa_order)

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

    return {
        "records":records,
        "events":events,
        "final_basis":basis,
        "final_coefficients":C,
        "graph":manager.graph,
        "registry":manager.registry,
        "lineage":lineage,
        "final_overlap":S.copy(),
        "final_hamiltonian":H.copy(),
        "settings":{
            "dt":float(dt),
            "steps":int(steps),
            "spa_order":int(spa_order),
            "spawn_action_threshold":float(spawn_action_threshold),
            "overlap_block":float(overlap_block),
            "child_overlap_block":float(child_overlap_block),
            "max_basis":int(max_basis),
            "max_generation":int(max_generation),
            "children_per_event":int(children_per_event),
            "allow_repeated_spawning":bool(allow_repeated_spawning),
            "minimum_spawn_separation_steps":int(minimum_spawn_separation_steps),
            "position_shifts":[float(x) for x in position_shifts],
            "width_scales":[float(x) for x in width_scales],
            "momentum_directions":list(momentum_directions),
        },
    }
