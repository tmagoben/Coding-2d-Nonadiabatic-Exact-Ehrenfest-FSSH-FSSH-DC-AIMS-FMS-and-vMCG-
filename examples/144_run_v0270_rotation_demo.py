"""Demonstrate arbitrary-rotation covariance of a correlated TDVP step."""

import sys
from pathlib import Path

import numpy as np


root = Path(__file__).resolve().parents[1]
if str(root) not in sys.path:
    sys.path.insert(0, str(root))

from gaussian_dynamics import (
    CorrelatedGaussianSpinorStateV270,
    build_correlated_metric_system_v270,
    correlated_implicit_midpoint_step_v270,
    pack_correlated_parameters_v270,
    rotate_correlated_velocity_v270,
    two_state_ci_soc_model_v260,
)


if __name__ == "__main__":
    state = CorrelatedGaussianSpinorStateV270(
        q=[[-0.25, 0.10]], p=[[3.0, 0.20]],
        width_matrices=[[[1.70, 0.25], [0.25, 2.60]]],
        chirp_matrices=[[[0.02, 0.03], [0.03, -0.01]]],
        coefficients=[[1.0, 0.0]],
    ).normalized()
    model = two_state_ci_soc_model_v260()
    angle = 0.371
    rotation = np.asarray(
        [[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]]
    )
    transformed_state = state.coordinate_rotated(rotation)
    transformed_model = model.coordinate_rotated(rotation)
    base_system = build_correlated_metric_system_v270(state, model)
    transformed_system = build_correlated_metric_system_v270(
        transformed_state, transformed_model
    )
    velocity_error = np.max(
        np.abs(
            rotate_correlated_velocity_v270(state, base_system.velocity, rotation)
            - transformed_system.velocity
        )
    )
    base_step = correlated_implicit_midpoint_step_v270(state, model, 0.002)
    transformed_step = correlated_implicit_midpoint_step_v270(
        transformed_state, transformed_model, 0.002
    )
    step_error = np.max(
        np.abs(
            pack_correlated_parameters_v270(base_step.end.coordinate_rotated(rotation))
            - pack_correlated_parameters_v270(transformed_step.end)
        )
    )
    print(f"velocity covariance error: {velocity_error:.3e}")
    print(f"midpoint-step covariance error: {step_error:.3e}")
    print(f"width eigenvalues: {base_step.end.width_eigenvalues[0]}")
    assert velocity_error < 3.0e-8
    assert step_error < 3.0e-7
