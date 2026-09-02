"""Run matched exact-grid, one-packet, and controlled 2D CI+SOC trajectories."""

import os
import sys
from pathlib import Path


for _thread_key in (
    "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
    "BLIS_NUM_THREADS", "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
):
    os.environ[_thread_key] = "1"

root = Path(__file__).resolve().parents[1]
if str(root) not in sys.path:
    sys.path.insert(0, str(root))

import numpy as np

from gaussian_dynamics import (
    DiagonalGaussianSpinorStateV260,
    ExactGridSettingsV260,
    UniformGrid2DV260,
    multidimensional_state_on_grid_v260,
    normalize_spinor_grid_v260,
    phase_aligned_grid_error_v260,
    run_controlled_multidimensional_dynamics_v260,
    run_exact_grid_ci_soc_v260,
    run_multidimensional_tdvp_v260,
    two_state_ci_soc_model_v260,
)


if __name__ == "__main__":
    model = two_state_ci_soc_model_v260(
        mass_au=(50.0, 50.0),
        kappa=0.04,
        coupling=0.04,
        frequencies=(0.03, 0.03),
        soc_scale=0.01,
    )
    grid = UniformGrid2DV260.from_bounds((-6.0, 6.0), (-6.0, 6.0), (64, 64))
    initial = DiagonalGaussianSpinorStateV260(
        q=[[-0.25, 0.0]],
        p=[[3.0, 0.0]],
        widths=[[2.0, 2.0]],
        chirps=[[0.0, 0.0]],
        coefficients=[[1.0, 0.0]],
    ).normalized()
    psi0 = normalize_spinor_grid_v260(
        multidimensional_state_on_grid_v260(initial, grid), grid
    )
    exact = run_exact_grid_ci_soc_v260(
        model,
        grid,
        psi0,
        settings=ExactGridSettingsV260(dt_au=0.01, steps=6, store_every=6),
    )
    one = run_multidimensional_tdvp_v260(initial, model, 0.01, 6)
    controlled = run_controlled_multidimensional_dynamics_v260(
        initial, model, 0.01, 6
    )

    def error(state):
        candidate = normalize_spinor_grid_v260(
            multidimensional_state_on_grid_v260(state, grid), grid
        )
        return phase_aligned_grid_error_v260(exact.final_state, candidate, grid)

    print(f"one-packet error: {error(one.final_state):.9e}")
    print(f"controlled error: {error(controlled.final_state):.9e}")
    print(f"events: {controlled.event_counts}")
    print(f"maximum norm drift: {controlled.maximum_norm_drift:.3e}")
    assert np.isfinite(error(controlled.final_state))
