import numpy as np

from .dynamic_graph_aims import DynamicGraphTBF
from .managed_graph_aims import run_managed_graph_aims
from .exact_benchmark import run_exact_ci_reference
from .convergence import vector_refinement_study


def _initial_tbf(q0, p0, A, state=1):
    return DynamicGraphTBF(
        uid=0,
        state=int(state),
        q=np.asarray(q0, dtype=float),
        p=np.asarray(p0, dtype=float),
        A=np.asarray(A, dtype=float),
        node=("seed", 0),
    )


def run_managed_ci_case(
    dt,
    final_time,
    q0=np.array([0.55, 0.45]),
    p0=np.array([0.6, 0.8]),
    A=1.2*np.eye(2),
    state=1,
    spa_order=0,
    spawn_action_threshold=2e-6,
    max_basis=2,
):
    steps = int(round(final_time / dt))
    out = run_managed_graph_aims(
        [_initial_tbf(q0, p0, A, state=state)],
        [1.0 + 0.0j],
        dt=dt,
        steps=steps,
        spa_order=spa_order,
        spawn_action_threshold=spawn_action_threshold,
        overlap_block=0.9,
        max_basis=max_basis,
        condition_limit=1e10,
        store_every=steps,
    )
    return out


def compare_managed_to_exact(
    managed_dt=2e-4,
    exact_dt=5e-4,
    final_time=0.01,
    spa_order=0,
    grid_n=40,
):
    exact = run_exact_ci_reference(
        grid_n=grid_n,
        dt=exact_dt,
        final_time=final_time,
    )
    managed = run_managed_ci_case(
        managed_dt,
        final_time,
        spa_order=spa_order,
    )

    p_exact = exact["final_populations_adiabatic"]
    p_managed = managed["records"][-1]["state_populations"]

    return {
        "exact_populations": p_exact,
        "managed_populations": p_managed,
        "population_l2_error": float(np.linalg.norm(p_managed - p_exact)),
        "exact_norm": exact["norm"],
        "managed_norm": managed["records"][-1]["norm"],
        "spawn_events": [e for e in managed["events"] if e["kind"] == "spawn"],
    }


def managed_timestep_refinement(
    dts=(8e-4, 4e-4, 2e-4),
    final_time=0.008,
    spa_order=0,
):
    populations = []
    runs = []
    for dt in dts:
        run = run_managed_ci_case(
            dt,
            final_time,
            spa_order=spa_order,
        )
        runs.append(run)
        populations.append(run["records"][-1]["state_populations"])

    study = vector_refinement_study(np.asarray(dts), populations)
    return {
        "dts": np.asarray(dts),
        "populations": np.asarray(populations),
        "successive_errors": study.successive_errors,
        "observed_orders": study.observed_orders,
        "runs": runs,
    }


def spa_order_comparison(dt=2e-4, final_time=0.008):
    out0 = run_managed_ci_case(dt, final_time, spa_order=0)
    out1 = run_managed_ci_case(dt, final_time, spa_order=1)

    p0 = out0["records"][-1]["state_populations"]
    p1 = out1["records"][-1]["state_populations"]

    return {
        "spa0_populations": p0,
        "spa1_populations": p1,
        "difference_l2": float(np.linalg.norm(p1-p0)),
        "spa0": out0,
        "spa1": out1,
    }
