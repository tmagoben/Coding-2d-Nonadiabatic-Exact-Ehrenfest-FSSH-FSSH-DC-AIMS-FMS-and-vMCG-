from dataclasses import dataclass, field
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
    _child_momentum,
)
from .graph_gaussian import GraphGaussianTBF
from .gaussian_nd import analytic_overlap_equal_width
from .moving_graph_gaussian import (
    nuclear_seed_basis_time_matrix,
    metric_compatible_basis_connection,
    moving_basis_coefficient_step,
    generalized_norm,
)
from .spa_matrix_elements import build_graph_gaussian_matrices_spa, spa1_correction_norm


def _ensure_graph(step, basis, manager, provider):
    for b in basis:
        node = _center_node(step, b.uid)
        if node not in manager.frames:
            connect = [b.node] if b.node in manager.frames and b.node != node else []
            manager.add_from_provider(node, b.q, provider, connect_to=connect)
        b.node = node

    refs = {}
    for i in range(len(basis)):
        refs[(i, i)] = basis[i].node
        for j in range(i + 1, len(basis)):
            node = _centroid_node(step, basis[i].uid, basis[j].uid)
            if node not in manager.frames:
                manager.add_from_provider(
                    node,
                    0.5 * (basis[i].q + basis[j].q),
                    provider,
                    connect_to=[basis[i].node, basis[j].node],
                )
            refs[(i, j)] = node
            refs[(j, i)] = node
    return refs


def _graph_basis(basis, dimension):
    out = []
    for b in basis:
        coeff = np.zeros(dimension, dtype=complex)
        coeff[b.state] = 1.0
        out.append(GraphGaussianTBF(b.node, b.q.copy(), b.p.copy(), b.A.copy(), coeff))
    return out


def maybe_spawn_integrated(
    basis,
    provider,
    next_uid,
    tracker,
    dt,
    overlap_block=0.85,
    allow_repeated_spawning=False,
    last_spawn_step=None,
    current_step=0,
    minimum_spawn_separation_steps=0,
):
    """Timestep-aware spawn rule based on accumulated |v.d| dt exposure."""
    ns = provider.evaluate(basis[0].q).frame.shape[1]

    for parent in basis:
        for target in range(ns):
            if target == parent.state:
                continue
            if (not allow_repeated_spawning) and target in parent.spawned_targets:
                continue

            if allow_repeated_spawning and last_spawn_step is not None:
                previous_step = last_spawn_step.get((int(parent.uid), int(target)))
                if previous_step is not None:
                    if int(current_step) - int(previous_step) < int(minimum_spawn_separation_steps):
                        continue

            rate = _indicator(parent, target, provider)
            key = (int(parent.uid), int(target))
            ready, exposure = tracker.update(key, rate, dt)
            if not ready:
                continue

            p_child = _child_momentum(parent, target, provider)
            if p_child is None:
                continue

            redundant = False
            for existing in basis:
                if existing.state != target:
                    continue
                if not np.allclose(existing.A, parent.A, atol=1e-12):
                    continue
                overlap = abs(
                    analytic_overlap_equal_width(
                        existing.q,
                        existing.p,
                        parent.q,
                        p_child,
                        parent.A,
                    )
                )
                if overlap >= overlap_block:
                    redundant = True
                    break

            if redundant:
                continue

            if not allow_repeated_spawning:
                parent.spawned_targets.add(target)
            if last_spawn_step is not None:
                last_spawn_step[(int(parent.uid), int(target))] = int(current_step)
            tracker.consume(key)
            child = DynamicGraphTBF(
                next_uid,
                target,
                parent.q.copy(),
                p_child,
                parent.A.copy(),
                parent.node,
            )
            return child, parent.uid, float(exposure), float(rate)

    return None, None, None, None


def run_managed_graph_aims(
    initial_basis,
    C0,
    provider=None,
    dt=2.0e-4,
    steps=40,
    spa_order=0,
    spawn_action_threshold=2.0e-4,
    spawn_coupling_floor=1.0e-8,
    overlap_block=0.90,
    max_basis=6,
    condition_limit=1.0e8,
    eigenvalue_floor=1.0e-9,
    max_pruning_loss=1.0e-8,
    allow_repeated_spawning=False,
    minimum_spawn_separation_steps=0,
    store_every=5,
):
    """v0.9 graph-AIMS prototype with SPA order, adaptive spawning, and pruning."""
    if spa_order not in (0, 1):
        raise ValueError("spa_order must be 0 or 1")

    provider = provider or AnalyticCI2DFrameProvider()
    dimension = provider.evaluate(initial_basis[0].q).frame.shape[1]
    manager = IncrementalElectronicGraph(dimension)

    basis = [
        DynamicGraphTBF(
            b.uid,
            b.state,
            b.q.copy(),
            b.p.copy(),
            b.A.copy(),
            b.node,
            set(b.spawned_targets),
        )
        for b in initial_basis
    ]
    C = np.asarray(C0, dtype=complex).copy()
    if len(C) != len(basis):
        raise ValueError("C0 length mismatch")

    refs = _ensure_graph(0, basis, manager, provider)
    gbasis = _graph_basis(basis, dimension)
    mass_matrix = provider.evaluate(basis[0].q).mass_matrix

    def build_matrices(gb, refmap, order):
        return build_graph_gaussian_matrices_spa(
            gb,
            manager.registry,
            mass_matrix,
            reference_selector=lambda i, j: refmap[(i, j)],
            order=order,
        )

    S, H = build_matrices(gbasis, refs, spa_order)

    initial_events = []
    if len(basis) > 1 and np.linalg.cond(S) > condition_limit:
        pruning = prune_redundant_basis(
            C,
            S,
            condition_limit=condition_limit,
            eigenvalue_floor=eigenvalue_floor,
            max_projection_loss=max_pruning_loss,
        )
        if pruning.removed:
            removed_uids = [basis[i].uid for i in pruning.removed]
            basis = [basis[i] for i in pruning.keep]
            C = pruning.coefficients
            refs = _ensure_graph(0, basis, manager, provider)
            gbasis = _graph_basis(basis, dimension)
            S, H = build_matrices(gbasis, refs, spa_order)
            initial_events.append({
                "kind": "prune",
                "step": 0,
                "time": 0.0,
                "removed_uids": [int(u) for u in removed_uids],
                "projection_loss": float(pruning.projection_loss),
                "condition_before": float(pruning.condition_before),
                "condition_after": float(pruning.condition_after),
            })

    C /= np.sqrt(generalized_norm(C, S))

    exposure = CouplingExposureTracker(
        action_threshold=spawn_action_threshold,
        coupling_floor=spawn_coupling_floor,
    )
    next_uid = max(b.uid for b in basis) + 1
    last_spawn_step = {}
    events = list(initial_events)
    records = []

    def record(step):
        norm = generalized_norm(C, S)
        populations = np.zeros(dimension)
        for state in range(dimension):
            idx = [i for i, b in enumerate(basis) if b.state == state]
            if idx:
                block = S[np.ix_(idx, idx)]
                cc = C[idx]
                populations[state] = np.real(np.vdot(cc, block @ cc))
        if norm > 0.0:
            populations /= norm

        # SPA1 correction diagnostic can be evaluated without changing propagation.
        spa_delta = 0.0
        if len(basis):
            H0 = build_matrices(_graph_basis(basis, dimension), refs, 0)[1]
            H1 = build_matrices(_graph_basis(basis, dimension), refs, 1)[1]
            spa_delta = spa1_correction_norm(H0, H1)

        records.append({
            "step": int(step),
            "time": float(step * dt),
            "norm": float(norm),
            "basis_size": int(len(basis)),
            "state_populations": populations.copy(),
            "condition_number": float(np.linalg.cond(S)),
            "spa_order": int(spa_order),
            "spa1_relative_correction": float(spa_delta),
            "graph": manager.summary(),
        })

    record(0)

    for step in range(1, steps + 1):
        old = [
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
        old_refs = refs
        old_gbasis = _graph_basis(old, dimension)
        S0, H0 = S, H

        for b in basis:
            b.q, b.p = _verlet(b, provider, dt)

        refs = _ensure_graph(step, basis, manager, provider)
        gbasis = _graph_basis(basis, dimension)
        S1, H1 = build_matrices(gbasis, refs, spa_order)

        qdots = []
        pdots = []
        for b in old:
            qdot, pdot = _kinematics(b, provider)
            qdots.append(qdot)
            pdots.append(pdot)

        seed = nuclear_seed_basis_time_matrix(
            old_gbasis,
            manager.registry,
            lambda i, j: old_refs[(i, j)],
            np.asarray(qdots),
            np.asarray(pdots),
        )
        T = metric_compatible_basis_connection(S0, S1, dt, seed=seed)
        C = moving_basis_coefficient_step(C, S0, H0, S1, H1, T, dt)
        S, H = S1, H1

        if len(basis) < max_basis:
            child, parent_uid, action, instantaneous_rate = maybe_spawn_integrated(
                basis,
                provider,
                next_uid,
                exposure,
                dt,
                overlap_block=overlap_block,
                allow_repeated_spawning=allow_repeated_spawning,
                last_spawn_step=last_spawn_step,
                current_step=step,
                minimum_spawn_separation_steps=minimum_spawn_separation_steps,
            )
            if child is not None:
                basis.append(child)
                C = np.concatenate([C, [0.0 + 0.0j]])
                events.append({
                    "kind": "spawn",
                    "step": int(step),
                    "time": float(step * dt),
                    "parent_uid": int(parent_uid),
                    "child_uid": int(next_uid),
                    "target_state": int(child.state),
                    "integrated_coupling_action": float(action),
                    "instantaneous_coupling_rate": float(instantaneous_rate),
                })
                next_uid += 1

                refs = _ensure_graph(step, basis, manager, provider)
                gbasis = _graph_basis(basis, dimension)
                S, H = build_matrices(gbasis, refs, spa_order)

        if len(basis) > 1 and np.linalg.cond(S) > condition_limit:
            pruning = prune_redundant_basis(
                C,
                S,
                condition_limit=condition_limit,
                eigenvalue_floor=eigenvalue_floor,
                max_projection_loss=max_pruning_loss,
            )

            if pruning.removed:
                removed_uids = [basis[i].uid for i in pruning.removed]
                basis = [basis[i] for i in pruning.keep]
                C = pruning.coefficients
                refs = _ensure_graph(step, basis, manager, provider)
                gbasis = _graph_basis(basis, dimension)
                S, H = build_matrices(gbasis, refs, spa_order)

                # Projection may reduce norm slightly; retain the projected physical
                # state rather than silently renormalizing it.  The loss is reported.
                events.append({
                    "kind": "prune",
                    "step": int(step),
                    "time": float(step * dt),
                    "removed_uids": [int(u) for u in removed_uids],
                    "projection_loss": float(pruning.projection_loss),
                    "condition_before": float(pruning.condition_before),
                    "condition_after": float(pruning.condition_after),
                })

        if step % store_every == 0:
            record(step)

    return {
        "records": records,
        "events": events,
        "final_basis": basis,
        "final_coefficients": C,
        "graph": manager.graph,
        "registry": manager.registry,
        "settings": {
            "dt": float(dt),
            "steps": int(steps),
            "spa_order": int(spa_order),
            "spawn_action_threshold": float(spawn_action_threshold),
            "condition_limit": float(condition_limit),
            "allow_repeated_spawning": bool(allow_repeated_spawning),
            "minimum_spawn_separation_steps": int(minimum_spawn_separation_steps),
        },
    }
