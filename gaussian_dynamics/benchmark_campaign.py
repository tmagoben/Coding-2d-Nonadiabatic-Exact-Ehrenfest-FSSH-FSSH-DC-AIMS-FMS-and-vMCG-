from dataclasses import dataclass, asdict
from itertools import product
import numpy as np

from .dynamic_graph_aims import DynamicGraphTBF
from .dynamic_gauge_graph import AnalyticCI2DFrameProvider
from .managed_graph_aims import run_managed_graph_aims
from .exact_benchmark import run_exact_ci_reference
from .benchmark_metrics import summarize_managed_run


@dataclass(frozen=True)
class CIPassageConfig:
    """A finite-impact-parameter passage near the 2D conical intersection."""
    q0: tuple = (-0.60, 0.25)
    p0: tuple = (10.0, 0.0)
    A_diag: tuple = (1.4, 1.4)
    state: int = 1
    mass: float = 5.0
    final_time: float = 0.60
    half_width: float = 4.0

    def q_array(self):
        return np.asarray(self.q0, dtype=float)

    def p_array(self):
        return np.asarray(self.p0, dtype=float)

    def A_matrix(self):
        return np.diag(np.asarray(self.A_diag, dtype=float))


def _initial_tbf(config):
    return DynamicGraphTBF(
        uid=0,
        state=int(config.state),
        q=config.q_array(),
        p=config.p_array(),
        A=config.A_matrix(),
        node=("seed", 0),
    )


def run_managed_passage(
    config=CIPassageConfig(),
    dt=0.005,
    spa_order=0,
    spawn_action_threshold=2e-4,
    max_basis=4,
    overlap_block=0.90,
    condition_limit=1e8,
    eigenvalue_floor=1e-9,
    max_pruning_loss=1e-8,
    allow_repeated_spawning=True,
    minimum_spawn_separation_steps=8,
    store_every=None,
):
    steps = int(round(config.final_time / dt))
    if steps <= 0:
        raise ValueError("final_time/dt must give at least one step.")
    if store_every is None:
        store_every = max(1, steps // 20)

    provider = AnalyticCI2DFrameProvider(nuclear_mass_au=config.mass)

    return run_managed_graph_aims(
        [_initial_tbf(config)],
        [1.0 + 0.0j],
        provider=provider,
        dt=dt,
        steps=steps,
        spa_order=spa_order,
        spawn_action_threshold=spawn_action_threshold,
        overlap_block=overlap_block,
        max_basis=max_basis,
        condition_limit=condition_limit,
        eigenvalue_floor=eigenvalue_floor,
        max_pruning_loss=max_pruning_loss,
        allow_repeated_spawning=allow_repeated_spawning,
        minimum_spawn_separation_steps=minimum_spawn_separation_steps,
        store_every=store_every,
    )


def run_exact_passage(
    config=CIPassageConfig(),
    grid_n=64,
    dt=0.0025,
):
    return run_exact_ci_reference(
        q0=config.q_array(),
        p0=config.p_array(),
        A=config.A_matrix(),
        state=config.state,
        mass=config.mass,
        grid_n=int(grid_n),
        half_width=config.half_width,
        dt=float(dt),
        final_time=config.final_time,
    )


def run_managed_parameter_surface(
    config=CIPassageConfig(),
    dts=(0.01, 0.005),
    spa_orders=(0, 1),
    spawn_action_thresholds=(1e-4, 2e-4),
    max_basis_values=(2, 4),
    overlap_blocks=(0.90,),
    exact_reference_populations=None,
):
    """Cartesian-product convergence/sensitivity surface for managed dynamics."""
    rows = []

    for dt, spa, spawn, max_basis, overlap_block in product(
        dts,
        spa_orders,
        spawn_action_thresholds,
        max_basis_values,
        overlap_blocks,
    ):
        run = run_managed_passage(
            config=config,
            dt=float(dt),
            spa_order=int(spa),
            spawn_action_threshold=float(spawn),
            max_basis=int(max_basis),
            overlap_block=float(overlap_block),
        )

        metrics = summarize_managed_run(run, exact_reference_populations)

        rows.append({
            "dt": float(dt),
            "spa_order": int(spa),
            "spawn_action_threshold": float(spawn),
            "max_basis": int(max_basis),
            "overlap_block": float(overlap_block),
            **metrics.to_dict(),
        })

    return rows


def run_exact_grid_timestep_surface(
    config=CIPassageConfig(),
    grid_values=(40, 56, 72),
    dt_values=(0.01, 0.005, 0.0025),
):
    """Exact-grid discretization surface.

    Each result is a fully independent split-operator calculation.  The finest
    `(grid_n, dt)` pair can be used as the numerical reference after checking the
    neighboring points rather than assuming it is converged.
    """
    rows = []

    for grid_n, dt in product(grid_values, dt_values):
        out = run_exact_passage(config=config, grid_n=int(grid_n), dt=float(dt))
        rows.append({
            "grid_n": int(grid_n),
            "dt": float(dt),
            "norm": float(out["norm"]),
            "populations": np.asarray(
                out["final_populations_adiabatic"], dtype=float
            ),
        })

    return rows


def select_finest_exact_reference(rows):
    """Select smallest dt and largest grid from a completed exact surface."""
    if not rows:
        raise ValueError("rows cannot be empty.")

    best = sorted(rows, key=lambda r: (float(r["dt"]), -int(r["grid_n"])))[0]
    return best


def campaign_settings_dict(config):
    return asdict(config)
