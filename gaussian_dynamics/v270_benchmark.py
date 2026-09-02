"""Cumulative acceptance campaign for the v0.27.0 correlated-width release."""

from dataclasses import asdict, dataclass, replace

import numpy as np

from .correlated_basis_adaptation_v270 import (
    ControlledCorrelatedBasisSettingsV270,
    V270_MULTIDIMENSIONAL_BASIS_CLAIMS,
    adapt_correlated_basis_once_v270,
    generate_correlated_spawn_candidates_v270,
    project_correlated_state_v270,
    run_controlled_correlated_dynamics_v270,
)
from .correlated_gaussian_tdvp_v270 import (
    CorrelatedGaussianSpinorStateV270,
    CorrelatedVariationalSettingsV270,
    V270_CORRELATED_TDVP_CLAIMS,
    build_correlated_metric_system_v270,
    correlated_implicit_midpoint_step_v270,
    smat_v270,
    svec_v270,
)
from .correlated_validation_v270 import run_correlated_validation_evidence_v270
from .multidimensional_soc_v260 import two_state_ci_soc_model_v260
from .v260_benchmark import run_v0260_release_benchmark


@dataclass(frozen=True)
class V270AcceptanceThresholds:
    expected_inherited_gates: int = 825
    expected_validation_gates: int = 100
    expected_core_gates: int = 35
    expected_new_gates: int = 135
    expected_total_gates: int = 960


def _raises_v270(callable_object, exceptions, text=None):
    try:
        callable_object()
    except exceptions as exc:
        return text is None or text in str(exc)
    return False


def _state_v270():
    return CorrelatedGaussianSpinorStateV270(
        [[-0.25, 0.10]], [[3.0, 0.20]],
        [[[1.70, 0.25], [0.25, 2.60]]],
        [[[0.02, 0.03], [0.03, -0.01]]],
        [[1.0, 0.0]],
    ).normalized()


def _core_controls_v270(evidence):
    settings = ControlledCorrelatedBasisSettingsV270().validate()
    tdvp_settings = settings.tdvp_settings
    state = _state_v270()
    model = two_state_ci_soc_model_v260()
    spawn_event = adapt_correlated_basis_once_v270(state, model)
    metric = build_correlated_metric_system_v270(state, model)
    step = correlated_implicit_midpoint_step_v270(state, model, 0.002)
    controlled = run_controlled_correlated_dynamics_v270(state, model, 0.001, 2)
    capped = adapt_correlated_basis_once_v270(
        state, model, settings=replace(settings, maximum_packet_count=1)
    )
    duplicate_q = np.vstack((state.q, state.q))
    duplicate_p = np.vstack((state.p, state.p))
    duplicate_width = np.concatenate((state.width_matrices, state.width_matrices), axis=0)
    duplicate_chirp = np.concatenate((state.chirp_matrices, state.chirp_matrices), axis=0)

    controls = {
        "disabled_full_width_contract_is_rejected": _raises_v270(
            lambda: replace(tdvp_settings, full_correlated_width_matrices=False).validate(),
            (ValueError,), "full-width variational contract",
        ),
        "disabled_orthogonal_covariance_contract_is_rejected": _raises_v270(
            lambda: replace(tdvp_settings, arbitrary_orthogonal_coordinate_covariance=False).validate(),
            (ValueError,), "full-width variational contract",
        ),
        "moving_electronic_frame_request_is_rejected": _raises_v270(
            lambda: replace(tdvp_settings, fixed_electronic_frame=False).validate(),
            (ValueError,), "full-width variational contract",
        ),
        "live_molecular_soc_request_is_rejected": _raises_v270(
            lambda: replace(tdvp_settings, real_molecular_soc_provider=True).validate(),
            (ValueError,), "does not admit live molecular-SOC trajectories",
        ),
        "invalid_width_domain_is_rejected": _raises_v270(
            lambda: replace(
                tdvp_settings, minimum_width_eigenvalue=2.0,
                maximum_width_eigenvalue=1.0,
            ).validate(),
            (ValueError,), "must exceed",
        ),
        "nonpositive_width_matrix_is_rejected": _raises_v270(
            lambda: CorrelatedGaussianSpinorStateV270(
                state.q, state.p, [[[1.0, 0.0], [0.0, -1.0]]],
                state.chirp_matrices, state.coefficients,
            ).validate(),
            (ValueError,), "positive definite",
        ),
        "nonsymmetric_width_matrix_is_rejected": _raises_v270(
            lambda: CorrelatedGaussianSpinorStateV270(
                state.q, state.p, [[[1.0, 0.3], [0.0, 1.0]]],
                state.chirp_matrices, state.coefficients,
            ).validate(),
            (ValueError,), "symmetric",
        ),
        "nonsymmetric_chirp_matrix_is_rejected": _raises_v270(
            lambda: CorrelatedGaussianSpinorStateV270(
                state.q, state.p, state.width_matrices,
                [[[0.0, 0.3], [0.0, 0.0]]], state.coefficients,
            ).validate(),
            (ValueError,), "symmetric",
        ),
        "invalid_matrix_array_shape_is_rejected": _raises_v270(
            lambda: CorrelatedGaussianSpinorStateV270(
                state.q, state.p, [[1.0, 2.0]], [[0.0, 0.0]], state.coefficients,
            ).validate(),
            (ValueError,), "shape",
        ),
        "invalid_active_shape_mask_is_rejected": _raises_v270(
            lambda: build_correlated_metric_system_v270(
                state, model, active_shape_mask=[True, False]
            ),
            (ValueError,), "one value per packet",
        ),
        "projection_vector_shaped_widths_are_rejected": _raises_v270(
            lambda: project_correlated_state_v270(
                state, state.q, state.p, [[1.0, 2.0]], [[0.0, 0.0]], model
            ),
            (ValueError,), "incompatible",
        ),
        "projection_nonsymmetric_width_is_rejected": _raises_v270(
            lambda: project_correlated_state_v270(
                state, state.q, state.p, [[[1.0, 0.3], [0.0, 1.0]]],
                state.chirp_matrices, model,
            ),
            (ValueError,), "symmetric",
        ),
        "projection_nonpositive_width_is_rejected": _raises_v270(
            lambda: project_correlated_state_v270(
                state, state.q, state.p, [[[1.0, 0.0], [0.0, -1.0]]],
                state.chirp_matrices, model,
            ),
            (ValueError,), "positive definite",
        ),
        "duplicate_projection_basis_is_rejected": _raises_v270(
            lambda: project_correlated_state_v270(
                state, duplicate_q, duplicate_p, duplicate_width, duplicate_chirp, model
            ),
            (ValueError,), "rank deficient",
        ),
        "isotropic_principal_axes_fail_closed": len(
            generate_correlated_spawn_candidates_v270(
                CorrelatedGaussianSpinorStateV270(
                    state.q, state.p, [2.0 * np.eye(2)], state.chirp_matrices,
                    state.coefficients,
                ).normalized()
            )
        ) == 0,
        "tampered_spawn_direction_policy_is_rejected": _raises_v270(
            lambda: replace(settings, spawn_directions="laboratory axes").validate(),
            (ValueError,), "spawn directions",
        ),
        "invalid_principal_axis_gap_is_rejected": _raises_v270(
            lambda: replace(settings, minimum_principal_axis_relative_gap=1.0).validate(),
            (ValueError,), "below one",
        ),
        "multiple_events_request_is_rejected": _raises_v270(
            lambda: replace(settings, maximum_events_per_checkpoint=2).validate(),
            (ValueError,), "at most one event",
        ),
        "disabled_lifecycle_procedure_is_rejected": _raises_v270(
            lambda: replace(settings, spawn_enabled=False).validate(),
            (ValueError,), "freezes spawn, prune, and merge",
        ),
        "maximum_packet_count_fails_closed": (
            capped.event_kind == "none" and capped.reason == "maximum packet count reached"
        ),
        "packet_serial_collision_is_rejected": _raises_v270(
            lambda: adapt_correlated_basis_once_v270(
                state, model, packet_ids=("g000001",), packet_ages=(0,),
                next_packet_serial=1,
            ),
            (ValueError,), "collides",
        ),
        "negative_packet_age_is_rejected": _raises_v270(
            lambda: adapt_correlated_basis_once_v270(
                state, model, packet_ids=("a",), packet_ages=(-1,), next_packet_serial=1
            ),
            (ValueError,), "packet age",
        ),
        "duplicate_packet_identity_is_rejected": _raises_v270(
            lambda: adapt_correlated_basis_once_v270(
                spawn_event.after, model, packet_ids=("a", "a"),
                packet_ages=(0, 0), next_packet_serial=2,
            ),
            (ValueError,), "unique",
        ),
        "nonlinear_nonconvergence_is_rejected": _raises_v270(
            lambda: correlated_implicit_midpoint_step_v270(
                state, model, 0.5,
                settings=replace(
                    tdvp_settings, nonlinear_max_function_evaluations=1
                ).validate(),
            ),
            (RuntimeError,), "implicit midpoint TDVP solve failed",
        ),
        "tampered_metric_velocity_is_rejected": _raises_v270(
            lambda: replace(metric, velocity=metric.velocity + 1.0e-3).validate(),
            (ValueError,), "disagrees with its SVD solve",
        ),
        "tampered_midpoint_residual_is_rejected": _raises_v270(
            lambda: replace(step, nonlinear_residual=step.nonlinear_residual + 1.0e-3).validate(),
            (ValueError,), "nonlinear residual",
        ),
        "tampered_projection_loss_is_rejected": _raises_v270(
            lambda: replace(
                spawn_event.projection,
                relative_projection_loss=spawn_event.projection.relative_projection_loss + 1.0e-3,
            ).validate(),
            (ValueError,), "relative projection loss",
        ),
        "tampered_spawn_identity_is_rejected": _raises_v270(
            lambda: replace(spawn_event, packet_ids_after=("g000000", "bad")).validate(),
            (ValueError,), "packet-ID order",
        ),
        "tampered_trajectory_metadata_is_rejected": _raises_v270(
            lambda: replace(
                controlled,
                final_packet_ids=tuple(
                    f"bad{index}" for index in range(controlled.final_state.ngaussian)
                ),
            ).validate(),
            (ValueError,), "final packet metadata",
        ),
        "inactive_matrix_shapes_are_bitwise_frozen": (
            lambda dormant, frozen_step: all(
                np.array_equal(left, right)
                for left, right in (
                    (frozen_step.start.q[-1], frozen_step.end.q[-1]),
                    (frozen_step.start.p[-1], frozen_step.end.p[-1]),
                    (frozen_step.start.width_matrices[-1], frozen_step.end.width_matrices[-1]),
                    (frozen_step.start.chirp_matrices[-1], frozen_step.end.chirp_matrices[-1]),
                )
            )(
                spawn_event.after,
                correlated_implicit_midpoint_step_v270(
                    spawn_event.after, model, 0.001, active_shape_mask=[True, False]
                ),
            )
        ),
        "offdiagonal_state_cannot_downgrade_to_v0260": _raises_v270(
            lambda: state.to_diagonal_v260(),
            (ValueError,), "cannot be represented",
        ),
        "nonorthogonal_coordinate_transform_is_rejected": _raises_v270(
            lambda: state.coordinate_rotated([[1.0, 0.2], [0.0, 1.0]]),
            (ValueError,), "must be orthogonal",
        ),
        "svec_rejects_nonsymmetric_input": _raises_v270(
            lambda: svec_v270([[1.0, 0.3], [0.0, 1.0]]),
            (ValueError,), "symmetric",
        ),
        "smat_rejects_wrong_vector_length": _raises_v270(
            lambda: smat_v270([1.0, 2.0], 2),
            (ValueError,), "invalid shape",
        ),
        "canonical_validation_fingerprint_is_stable": bool(
            len(evidence.fingerprint()) == 64
            and evidence.fingerprint() == evidence.fingerprint()
            and all(character in "0123456789abcdef" for character in evidence.fingerprint())
        ),
    }
    if len(controls) != 35:
        raise AssertionError(f"v0.27.0 must define exactly 35 core gates, found {len(controls)}.")
    return {name: bool(value) for name, value in controls.items()}


def run_v0270_release_benchmark(
    thresholds=V270AcceptanceThresholds(), *, memory_probe_policy="proc_self"
):
    inherited = run_v0260_release_benchmark(memory_probe_policy=memory_probe_policy)
    evidence = run_correlated_validation_evidence_v270()
    if evidence.check_count != thresholds.expected_validation_gates:
        raise AssertionError("v0.27.0 validation evidence must define 100 gates.")
    core = _core_controls_v270(evidence)
    inherited_checks = {
        f"inherited_v0260::{name}": bool(value)
        for name, value in inherited["acceptance"]["checks"].items()
    }
    validation_checks = {
        f"correlated_validation::{name}": bool(value)
        for name, value in evidence.checks.items()
    }
    core_checks = {
        f"correlated_core::{name}": bool(value) for name, value in core.items()
    }
    if len(inherited_checks) != thresholds.expected_inherited_gates:
        raise AssertionError("v0.27.0 must inherit exactly 825 v0.26.0 gates.")
    new_checks = {**validation_checks, **core_checks}
    if len(new_checks) != thresholds.expected_new_gates:
        raise AssertionError("v0.27.0 must define exactly 135 new gates.")
    checks = {**inherited_checks, **new_checks}
    if len(checks) != thresholds.expected_total_gates:
        raise AssertionError("v0.27.0 must define exactly 960 cumulative gates.")
    return {
        "release": "v0.27.0",
        "theme": (
            "full complex-symmetric correlated Gaussian widths and chirps, exact "
            "multivariate moments, rotation-covariant McLachlan TDVP, and controlled "
            "intrinsic-axis basis adaptation"
        ),
        "correlated_validation_evidence": evidence.as_dict(),
        "correlated_core_controls": core,
        "claims": {
            "tdvp": dict(V270_CORRELATED_TDVP_CLAIMS),
            "basis": dict(V270_MULTIDIMENSIONAL_BASIS_CLAIMS),
        },
        "inherited_v0260": inherited,
        "acceptance": {
            "passed": bool(all(checks.values())),
            "checks": checks,
            "inherited_gate_count": len(inherited_checks),
            "validation_gate_count": len(validation_checks),
            "core_gate_count": len(core_checks),
            "new_gate_count": len(new_checks),
            "total_gate_count": len(checks),
            "thresholds": asdict(thresholds),
        },
    }


def save_v0270_release_benchmark(path, **kwargs):
    from .campaign_io import save_campaign_json

    result = run_v0270_release_benchmark(**kwargs)
    return save_campaign_json(path, result)
