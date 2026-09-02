"""Cumulative reference-first multidimensional release campaign for v0.26.0."""

from dataclasses import asdict, dataclass, replace

import numpy as np

from .multidimensional_basis_adaptation_v260 import (
    ControlledMultidimensionalBasisSettingsV260,
    MultidimensionalSpawnCandidateV260,
    V260_MULTIDIMENSIONAL_BASIS_CLAIMS,
    adapt_multidimensional_basis_once_v260,
    evaluate_multidimensional_spawn_candidate_v260,
    project_multidimensional_state_v260,
    run_controlled_multidimensional_dynamics_v260,
)
from .multidimensional_gaussian_tdvp_v260 import (
    DiagonalGaussianSpinorStateV260,
    V260_MULTIDIMENSIONAL_TDVP_CLAIMS,
    build_multidimensional_metric_system_v260,
    multidimensional_implicit_midpoint_step_v260,
)
from .multidimensional_soc_v260 import (
    ExactGridSettingsV260,
    QuadraticSpinHamiltonianNDV260,
    UniformGrid2DV260,
    V260_EXACT_GRID_CLAIMS,
    normalize_spinor_grid_v260,
    two_state_ci_soc_model_v260,
)
from .multidimensional_validation_v260 import (
    run_multidimensional_validation_evidence_v260,
)
from .v253_benchmark import run_v0253_release_benchmark


@dataclass(frozen=True)
class V260AcceptanceThresholds:
    expected_inherited_gates: int = 715
    expected_validation_gates: int = 80
    expected_core_gates: int = 30
    expected_new_gates: int = 110
    expected_total_gates: int = 825


def _raises_v260(callable_object, exceptions, text=None):
    try:
        callable_object()
    except exceptions as exc:
        return text is None or text in str(exc)
    return False


def _core_controls_v260(evidence):
    settings = ControlledMultidimensionalBasisSettingsV260().validate()
    model = two_state_ci_soc_model_v260()
    state = DiagonalGaussianSpinorStateV260(
        q=[[-2.0, 0.0]],
        p=[[7.0, 0.0]],
        widths=[[1.0, 1.0]],
        chirps=[[0.0, 0.0]],
        coefficients=[[1.0, 0.0]],
    ).normalized()
    spawn_event = adapt_multidimensional_basis_once_v260(state, model)
    metric = build_multidimensional_metric_system_v260(state, model)
    duplicate = MultidimensionalSpawnCandidateV260(
        state.q[0], state.p[0], state.widths[0], state.chirps[0],
        0, "position", 0, 1,
    )
    capped = adapt_multidimensional_basis_once_v260(
        state, model, settings=replace(settings, maximum_packet_count=1)
    )
    short_trajectory = run_controlled_multidimensional_dynamics_v260(
        state, model, 0.01, 2
    )
    controls = {
        "full_width_matrix_request_is_rejected": _raises_v260(
            lambda: replace(
                settings.tdvp_settings, full_correlated_width_matrices=True
            ).validate(),
            (ValueError,), "does not admit full correlated width matrices",
        ),
        "moving_electronic_frame_request_is_rejected": _raises_v260(
            lambda: replace(settings.tdvp_settings, fixed_electronic_frame=False).validate(),
            (ValueError,), "fixed electronic frame",
        ),
        "real_molecular_soc_request_is_rejected": _raises_v260(
            lambda: replace(settings.tdvp_settings, real_molecular_soc_provider=True).validate(),
            (ValueError,), "does not admit live molecular-SOC trajectories",
        ),
        "disabled_lifecycle_procedure_is_rejected": _raises_v260(
            lambda: replace(settings, spawn_enabled=False).validate(),
            (ValueError,), "freezes spawn, prune, and merge",
        ),
        "multiple_events_request_is_rejected": _raises_v260(
            lambda: replace(settings, maximum_events_per_checkpoint=2).validate(),
            (ValueError,), "at most one event",
        ),
        "invalid_novelty_gate_is_rejected": _raises_v260(
            lambda: replace(settings, minimum_candidate_novelty=1.0).validate(),
            (ValueError,), "below one",
        ),
        "invalid_merge_overlap_is_rejected": _raises_v260(
            lambda: replace(settings, minimum_merge_overlap=1.1).validate(),
            (ValueError,), "cannot exceed one",
        ),
        "nonpositive_mass_is_rejected": _raises_v260(
            lambda: QuadraticSpinHamiltonianNDV260(
                [[1.0, 0.0], [0.0, -1.0]], model.H0, model.H1, model.H2
            ).validate(),
            (ValueError,), "positive definite",
        ),
        "nonhermitian_hamiltonian_is_rejected": _raises_v260(
            lambda: QuadraticSpinHamiltonianNDV260(
                model.mass_matrix_au,
                model.H0 + np.asarray([[0.0, 0.1], [0.0, 0.0]]),
                model.H1,
                model.H2,
            ).validate(),
            (ValueError,), "must be Hermitian",
        ),
        "incomplete_projector_resolution_is_rejected": _raises_v260(
            lambda: QuadraticSpinHamiltonianNDV260(
                model.mass_matrix_au,
                model.H0,
                model.H1,
                model.H2,
                projectors={"partial": np.diag([1.0, 0.0])},
            ).validate(),
            (ValueError,), "resolve the identity",
        ),
        "nonorthogonal_coordinate_transform_is_rejected": _raises_v260(
            lambda: model.coordinate_rotated([[1.0, 0.2], [0.0, 1.0]]),
            (ValueError,), "must be orthogonal",
        ),
        "anisotropic_general_width_rotation_is_rejected": _raises_v260(
            lambda: DiagonalGaussianSpinorStateV260(
                [[0.0, 0.0]], [[0.0, 0.0]], [[1.0, 2.0]], [[0.0, 0.0]], [[1.0, 0.0]]
            ).normalized().coordinate_rotated(
                [[2.0**-0.5, -2.0**-0.5], [2.0**-0.5, 2.0**-0.5]]
            ),
            (ValueError,), "not closed under general rotations",
        ),
        "nonuniform_grid_is_rejected": _raises_v260(
            lambda: UniformGrid2DV260(
                [0.0, 1.0, 2.1, 3.0, 4.0, 5.0, 6.0, 7.0],
                np.arange(8.0),
            ).validate(),
            (ValueError,), "uniformly increasing",
        ),
        "undersized_grid_is_rejected": _raises_v260(
            lambda: UniformGrid2DV260.from_bounds((-1.0, 1.0), (-1.0, 1.0), (7, 8)),
            (ValueError,), "at least eight",
        ),
        "zero_grid_timestep_is_rejected": _raises_v260(
            lambda: ExactGridSettingsV260(dt_au=0.0).validate(),
            (ValueError,), "nonzero",
        ),
        "invalid_spinor_grid_shape_is_rejected": _raises_v260(
            lambda: normalize_spinor_grid_v260(
                np.zeros((2, 8, 7)),
                UniformGrid2DV260.from_bounds((-1.0, 1.0), (-1.0, 1.0), (8, 8)),
            ),
            (ValueError,), "shape",
        ),
        "negative_gaussian_width_is_rejected": _raises_v260(
            lambda: DiagonalGaussianSpinorStateV260(
                [[0.0, 0.0]], [[0.0, 0.0]], [[1.0, -1.0]], [[0.0, 0.0]], [[1.0, 0.0]]
            ).validate(),
            (ValueError,), "must be positive",
        ),
        "invalid_active_shape_mask_is_rejected": _raises_v260(
            lambda: build_multidimensional_metric_system_v260(
                state, model, active_shape_mask=[True, False]
            ),
            (ValueError,), "one value per packet",
        ),
        "duplicate_projection_basis_is_rejected": _raises_v260(
            lambda: project_multidimensional_state_v260(
                state,
                np.vstack((state.q, state.q)),
                np.vstack((state.p, state.p)),
                np.vstack((state.widths, state.widths)),
                np.vstack((state.chirps, state.chirps)),
                model,
            ),
            (ValueError,), "rank deficient",
        ),
        "duplicate_spawn_candidate_is_rejected": (
            evaluate_multidimensional_spawn_candidate_v260(state, model, duplicate).admitted
            is False
        ),
        "negative_packet_age_is_rejected": _raises_v260(
            lambda: adapt_multidimensional_basis_once_v260(
                state, model, packet_ids=("a",), packet_ages=(-1,), next_packet_serial=1
            ),
            (ValueError,), "packet age",
        ),
        "duplicate_packet_id_is_rejected": _raises_v260(
            lambda: adapt_multidimensional_basis_once_v260(
                spawn_event.after,
                model,
                packet_ids=("a", "a"),
                packet_ages=(0, 0),
                next_packet_serial=2,
            ),
            (ValueError,), "unique",
        ),
        "packet_serial_collision_is_rejected": _raises_v260(
            lambda: adapt_multidimensional_basis_once_v260(
                state,
                model,
                packet_ids=("g000001",),
                packet_ages=(0,),
                next_packet_serial=1,
            ),
            (ValueError,), "collides",
        ),
        "maximum_packet_count_fails_closed": (
            capped.event_kind == "none" and capped.reason == "maximum packet count reached"
        ),
        "nonlinear_nonconvergence_is_rejected": _raises_v260(
            lambda: multidimensional_implicit_midpoint_step_v260(
                state,
                model,
                0.5,
                settings=replace(
                    settings.tdvp_settings,
                    nonlinear_max_function_evaluations=1,
                ).validate(),
            ),
            (RuntimeError,), "implicit midpoint TDVP solve failed",
        ),
        "tampered_projection_loss_is_rejected": _raises_v260(
            lambda: replace(
                spawn_event.projection,
                relative_projection_loss=spawn_event.projection.relative_projection_loss + 1.0e-3,
            ).validate(),
            (ValueError,), "relative projection loss",
        ),
        "tampered_metric_velocity_is_rejected": _raises_v260(
            lambda: replace(metric, velocity=metric.velocity + 1.0e-3).validate(),
            (ValueError,), "disagrees with its SVD solve",
        ),
        "tampered_spawn_event_id_is_rejected": _raises_v260(
            lambda: replace(
                spawn_event, packet_ids_after=("g000000", "bad")
            ).validate(),
            (ValueError,), "packet-ID order",
        ),
        "tampered_trajectory_metadata_is_rejected": _raises_v260(
            lambda: replace(
                short_trajectory,
                final_packet_ids=tuple(
                    f"bad{index}" for index in range(short_trajectory.final_state.ngaussian)
                ),
            ).validate(),
            (ValueError,), "final packet metadata",
        ),
        "canonical_validation_fingerprint_is_stable": bool(
            len(evidence.fingerprint()) == 64
            and evidence.fingerprint() == evidence.fingerprint()
            and all(character in "0123456789abcdef" for character in evidence.fingerprint())
        ),
    }
    if len(controls) != 30:
        raise AssertionError(
            f"v0.26.0 must define exactly 30 core gates, found {len(controls)}."
        )
    return {name: bool(value) for name, value in controls.items()}


def run_v0260_release_benchmark(
    thresholds=V260AcceptanceThresholds(), *, memory_probe_policy="proc_self"
):
    inherited = run_v0253_release_benchmark(memory_probe_policy=memory_probe_policy)
    evidence = run_multidimensional_validation_evidence_v260()
    if evidence.check_count != thresholds.expected_validation_gates:
        raise AssertionError("v0.26.0 validation evidence must define 80 gates.")
    core = _core_controls_v260(evidence)
    inherited_checks = {
        f"inherited_v0253::{name}": bool(value)
        for name, value in inherited["acceptance"]["checks"].items()
    }
    validation_checks = {
        f"multidimensional_validation::{name}": bool(value)
        for name, value in evidence.checks.items()
    }
    core_checks = {
        f"multidimensional_core::{name}": bool(value)
        for name, value in core.items()
    }
    if len(inherited_checks) != thresholds.expected_inherited_gates:
        raise AssertionError("v0.26.0 must inherit exactly 715 v0.25.3 gates.")
    new_checks = {**validation_checks, **core_checks}
    if len(new_checks) != thresholds.expected_new_gates:
        raise AssertionError("v0.26.0 must define exactly 110 new gates.")
    checks = {**inherited_checks, **new_checks}
    if len(checks) != thresholds.expected_total_gates:
        raise AssertionError("v0.26.0 must define exactly 825 cumulative gates.")
    return {
        "release": "v0.26.0",
        "theme": (
            "reference-first two-dimensional CI+SOC dynamics, complete doublet/triplet "
            "models, diagonal-width multidimensional TDVP, and controlled basis adaptation"
        ),
        "multidimensional_validation_evidence": evidence.as_dict(),
        "multidimensional_core_controls": core,
        "claims": {
            "exact_grid": dict(V260_EXACT_GRID_CLAIMS),
            "tdvp": dict(V260_MULTIDIMENSIONAL_TDVP_CLAIMS),
            "basis": dict(V260_MULTIDIMENSIONAL_BASIS_CLAIMS),
        },
        "inherited_v0253": inherited,
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


def save_v0260_release_benchmark(path, **kwargs):
    from .campaign_io import save_campaign_json

    result = run_v0260_release_benchmark(**kwargs)
    return save_campaign_json(path, result)
