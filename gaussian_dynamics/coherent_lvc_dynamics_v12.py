import numpy as np

from .adaptive_spawning import CouplingExposureTracker
from .basis_management import prune_redundant_basis
from .dynamic_gauge_graph import AnalyticCI2DFrameProvider
from .dynamic_graph_aims import DynamicGraphTBF, _verlet, _kinematics
from .local_diabatic_tbf_v12 import (
    LocalDiabaticTBF,
    from_adiabatic_guided_tbf,
    reset_to_instantaneous_adiabatic_spinor,
    parallel_transport_spinor_full_space,
)
from .managed_graph_aims_v11 import _spawn_from_integrated_exposure
from .moving_graph_gaussian import (
    metric_compatible_basis_connection,
    moving_basis_coefficient_step,
    generalized_norm,
)
from .moving_basis_v12 import moving_basis_midpoint_cayley_step
from .lvc_exact_gaussian import (
    build_exact_lvc_gaussian_matrices,
    exact_lvc_basis_time_matrix,
    exact_lvc_basis_time_matrix_with_spinor_derivatives,
    center_spinor_time_derivative,
)
from .electronic_observables import (
    reduced_electronic_density_analytic_ci_diabatic,
)


def _copy_basis(basis, provider=None):
    out=[]
    for b in basis:
        if isinstance(b, LocalDiabaticTBF):
            out.append(b.copy())
        elif hasattr(b, "spinor"):
            out.append(
                LocalDiabaticTBF(
                    uid=b.uid,
                    state=b.state,
                    q=np.asarray(b.q,float).copy(),
                    p=np.asarray(b.p,float).copy(),
                    A=np.asarray(b.A,float).copy(),
                    spinor=np.asarray(b.spinor,complex).copy(),
                    node=getattr(b,"node",None),
                    spawned_targets=set(getattr(b,"spawned_targets",set())),
                )
            )
        else:
            if provider is None:
                raise ValueError("provider is required to initialize electronic spinors.")
            out.append(from_adiabatic_guided_tbf(b,provider))
    return out


def _midpoint_basis(old_basis, new_basis, provider, spinor_transport):
    if len(old_basis) != len(new_basis):
        raise ValueError("midpoint basis requires equal basis sizes.")

    out = []
    for old, new in zip(old_basis, new_basis):
        if old.uid != new.uid or old.state != new.state:
            raise ValueError("basis identity changed inside one propagation substep.")
        if not np.allclose(old.A, new.A, atol=1e-12):
            raise ValueError("v0.12 assumes each TBF width is frozen during a step.")

        qmid=0.5*(old.q+new.q)
        pmid=0.5*(old.p+new.p)

        if spinor_transport=="parallel":
            spinor=old.spinor.copy()
        elif spinor_transport=="instantaneous":
            point=provider.evaluate(qmid)
            spinor=np.asarray(point.frame[:,int(old.state)],complex)
        else:
            raise ValueError("unknown spinor transport.")

        out.append(
            LocalDiabaticTBF(
                uid=old.uid,
                state=old.state,
                q=qmid,
                p=pmid,
                A=old.A.copy(),
                spinor=spinor,
                node=old.node,
                spawned_targets=set(old.spawned_targets),
            )
        )
    return out

def _kinematic_arrays(basis, provider):
    qdots = []
    pdots = []
    for b in basis:
        qdot, pdot = _kinematics(b, provider)
        qdots.append(qdot)
        pdots.append(pdot)
    return np.asarray(qdots), np.asarray(pdots)


def _state_label_proxy(C, S, basis, nstate):
    """Legacy diagnostic only; not used as the rigorous electronic observable."""
    norm = generalized_norm(C, S)
    out = np.zeros(nstate, dtype=float)

    for state in range(nstate):
        idx = [i for i, b in enumerate(basis) if b.state == state]
        if idx:
            block = S[np.ix_(idx, idx)]
            cc = C[idx]
            out[state] = np.real(np.vdot(cc, block @ cc))

    if norm > 0.0:
        out /= norm
    return out


def run_coherent_lvc_gaussians(
    initial_basis,
    C0,
    provider=None,
    dt=0.005,
    steps=120,
    integrator="cayley",
    spinor_transport="parallel",
    spawn_action_threshold=1e-4,
    spawn_coupling_floor=1e-8,
    overlap_block=0.9999,
    child_overlap_block=0.995,
    max_basis=10,
    max_generation=5,
    children_per_event=2,
    allow_repeated_spawning=True,
    minimum_spawn_separation_steps=4,
    position_shifts=(0.0, 0.05, -0.05),
    width_scales=(0.65, 1.0, 1.55),
    momentum_directions=("nac", "momentum"),
    novelty_power=0.5,
    condition_limit=1e9,
    eigenvalue_floor=1e-10,
    max_pruning_loss=1e-7,
    store_every=20,
):
    r"""v0.12 coherence-oriented analytic-LVC Gaussian propagation.

    Nuclear centers
    ---------------
    TBF centers continue to move classically on their adiabatic surfaces.

    Electronic representation
    -------------------------
    Each Gaussian carries an explicit two-component spinor in the analytic model's
    global diabatic basis.  The default is overlap/parallel transport, which keeps the
    physical spinor locally diabatic while the nuclear center is still guided by an
    adiabatic surface.  `spinor_transport="instantaneous"` is retained as an ablation
    that resets the spinor to the adiabatic eigenvector at every center step.

    Hamiltonian
    -----------
    The full LVC potential matrix element is integrated analytically.  No SPA0/SPA1
    approximation is used in H for this benchmark-specialized propagator.

    Moving-basis connection
    -----------------------
    In the default parallel-transport mode the electronic spinor derivative
    vanishes in the complete global diabatic electronic space, so the seed contains
    the exact nuclear Gaussian derivative.  In the instantaneous-adiabatic ablation,
    the explicit electronic term from v . d is also included.

    The seed's anti-Hermitian content is retained while a minimal Hermitian correction
    enforces the discrete metric identity

        (S_{n+1}-S_n)/dt = T + T^dagger.

    This is a local-diabatic/frozen-spinor benchmark ansatz.  It is not claimed to be
    the full Born-Huang AIMS Hamiltonian used in production molecular calculations.
    """
    if integrator not in {"cayley", "rk4"}:
        raise ValueError("integrator must be 'cayley' or 'rk4'.")
    if spinor_transport not in {"parallel", "instantaneous"}:
        raise ValueError(
            "spinor_transport must be 'parallel' or 'instantaneous'."
        )
    if int(store_every) <= 0:
        raise ValueError("store_every must be positive.")

    provider = provider or AnalyticCI2DFrameProvider()
    if not isinstance(provider, AnalyticCI2DFrameProvider):
        raise TypeError(
            "v0.12 exact-LVC propagator is intentionally specialized to "
            "AnalyticCI2DFrameProvider."
        )

    basis = _copy_basis(initial_basis, provider=provider)
    C = np.asarray(C0, dtype=complex).copy()
    if len(C) != len(basis) or len(basis) == 0:
        raise ValueError("initial coefficient/basis sizes are inconsistent.")

    nstate = provider.evaluate(basis[0].q).frame.shape[1]

    lineage = {
        int(b.uid): {
            "parent_uid": None,
            "generation": 0,
            "birth_step": 0,
            "birth_time": 0.0,
            "birth_state": int(b.state),
            "width_scale_from_parent": 1.0,
        }
        for b in basis
    }

    S, H = build_exact_lvc_gaussian_matrices(basis, provider)

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
            removed = [basis[i].uid for i in pruning.removed]
            basis = [basis[i] for i in pruning.keep]
            C = pruning.coefficients
            S, H = build_exact_lvc_gaussian_matrices(basis, provider)
            initial_events.append({
                "kind": "prune",
                "step": 0,
                "time": 0.0,
                "removed_uids": [int(x) for x in removed],
                "projection_loss": float(pruning.projection_loss),
                "condition_before": float(pruning.condition_before),
                "condition_after": float(pruning.condition_after),
            })

    C /= np.sqrt(generalized_norm(C, S))

    tracker = CouplingExposureTracker(
        action_threshold=spawn_action_threshold,
        coupling_floor=spawn_coupling_floor,
    )
    next_uid = max(int(b.uid) for b in basis) + 1
    last_spawn_step = {}
    events = list(initial_events)
    records = []

    candidate_kwargs = dict(
        position_shifts=position_shifts,
        width_scales=width_scales,
        momentum_directions=momentum_directions,
        overlap_block=overlap_block,
        child_overlap_block=child_overlap_block,
        novelty_power=novelty_power,
    )

    latest_connection_correction = 0.0
    latest_connection_norm = 0.0

    def record(step):
        rho = reduced_electronic_density_analytic_ci_diabatic(
            C,
            basis,
            normalize=True,
        )
        generations = [lineage[int(b.uid)]["generation"] for b in basis]

        records.append({
            "step": int(step),
            "time": float(step*dt),
            "norm": float(generalized_norm(C, S)),
            "basis_size": int(len(basis)),
            "condition_number": float(np.linalg.cond(S)),
            "maximum_generation": int(max(generations, default=0)),
            "state_population_proxy": _state_label_proxy(C, S, basis, nstate),
            "diabatic_populations": np.real(np.diag(rho)),
            "diabatic_coherence": complex(rho[0, 1]),
            "connection_correction_norm": float(latest_connection_correction),
            "connection_seed_norm": float(latest_connection_norm),
        })

    record(0)

    for step in range(1, int(steps)+1):
        old_basis = _copy_basis(basis, provider=provider)
        S0 = S.copy()
        H0 = H.copy()

        qdot0, pdot0 = _kinematic_arrays(old_basis, provider)

        for old_b, b in zip(old_basis, basis):
            old_q=old_b.q.copy()
            b.q, b.p = _verlet(b, provider, dt)

            if spinor_transport=="parallel":
                parallel_transport_spinor_full_space(
                    b,
                    old_q,
                    b.q,
                    provider,
                )
            else:
                reset_to_instantaneous_adiabatic_spinor(
                    b,
                    provider,
                )

        qdot1, pdot1 = _kinematic_arrays(basis, provider)

        S1, H1 = build_exact_lvc_gaussian_matrices(basis, provider)

        middle = _midpoint_basis(
            old_basis,
            basis,
            provider,
            spinor_transport,
        )
        qdot_mid = 0.5*(qdot0+qdot1)
        pdot_mid = 0.5*(pdot0+pdot1)

        if spinor_transport=="parallel":
            T_seed = exact_lvc_basis_time_matrix(
                middle,
                provider,
                qdot_mid,
                pdot_mid,
            )
        else:
            spinor_dots=[
                center_spinor_time_derivative(
                    b,
                    provider,
                    qdot=qdot_mid[i],
                )
                for i,b in enumerate(middle)
            ]
            T_seed = exact_lvc_basis_time_matrix_with_spinor_derivatives(
                middle,
                provider,
                qdot_mid,
                pdot_mid,
                spinor_dots,
            )
        T_mid = metric_compatible_basis_connection(
            S0,
            S1,
            dt,
            seed=T_seed,
        )

        latest_connection_correction = np.linalg.norm(
            T_mid-T_seed, ord="fro"
        )
        latest_connection_norm = np.linalg.norm(T_seed, ord="fro")

        if integrator == "cayley":
            C = moving_basis_midpoint_cayley_step(
                C,
                S0,
                H0,
                S1,
                H1,
                T_mid,
                dt,
            )
        else:
            C = moving_basis_coefficient_step(
                C,
                S0,
                H0,
                S1,
                H1,
                T_mid,
                dt,
            )

        S, H = S1, H1

        if len(basis) < int(max_basis):
            spawn = _spawn_from_integrated_exposure(
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
                for child_raw, detail in zip(spawn["children"], spawn["details"]):
                    child = from_adiabatic_guided_tbf(
                        child_raw,
                        provider,
                    )
                    basis.append(child)
                    C = np.concatenate([C, [0.0+0.0j]])
                    lineage[int(child.uid)] = {
                        "parent_uid": int(spawn["parent_uid"]),
                        "generation": int(detail["generation"]),
                        "birth_step": int(step),
                        "birth_time": float(step*dt),
                        "birth_state": int(child.state),
                        "width_scale_from_parent": float(detail["width_scale"]),
                    }

                next_uid = int(spawn["next_uid"])
                events.append({
                    "kind": "optimized_spawn",
                    "step": int(step),
                    "time": float(step*dt),
                    "parent_uid": int(spawn["parent_uid"]),
                    "parent_generation": int(spawn["parent_generation"]),
                    "target_state": int(spawn["target_state"]),
                    "integrated_coupling_action": float(
                        spawn["integrated_coupling_action"]
                    ),
                    "instantaneous_coupling_rate": float(
                        spawn["instantaneous_coupling_rate"]
                    ),
                    "children": spawn["details"],
                })

                S, H = build_exact_lvc_gaussian_matrices(basis, provider)

        if len(basis) > 1 and np.linalg.cond(S) > condition_limit:
            pruning = prune_redundant_basis(
                C,
                S,
                condition_limit=condition_limit,
                eigenvalue_floor=eigenvalue_floor,
                max_projection_loss=max_pruning_loss,
            )

            if pruning.removed:
                removed = [basis[i].uid for i in pruning.removed]
                basis = [basis[i] for i in pruning.keep]
                C = pruning.coefficients
                S, H = build_exact_lvc_gaussian_matrices(basis, provider)

                events.append({
                    "kind": "prune",
                    "step": int(step),
                    "time": float(step*dt),
                    "removed_uids": [int(x) for x in removed],
                    "projection_loss": float(pruning.projection_loss),
                    "condition_before": float(pruning.condition_before),
                    "condition_after": float(pruning.condition_after),
                })

        if step % int(store_every) == 0:
            record(step)

    if records[-1]["step"] != int(steps):
        record(int(steps))

    return {
        "records": records,
        "events": events,
        "final_basis": basis,
        "final_coefficients": C,
        "final_overlap": S.copy(),
        "final_hamiltonian": H.copy(),
        "lineage": lineage,
        "settings": {
            "dt": float(dt),
            "steps": int(steps),
            "integrator": integrator,
            "spinor_transport": spinor_transport,
            "hamiltonian": "exact analytic 2D LVC frozen-spinor Gaussian matrix",
            "basis_connection": "nuclear + electronic spinor derivative",
            "spawn_action_threshold": float(spawn_action_threshold),
            "overlap_block": float(overlap_block),
            "child_overlap_block": float(child_overlap_block),
            "max_basis": int(max_basis),
            "max_generation": int(max_generation),
            "children_per_event": int(children_per_event),
            "position_shifts": [float(x) for x in position_shifts],
            "width_scales": [float(x) for x in width_scales],
            "momentum_directions": list(momentum_directions),
        },
    }
