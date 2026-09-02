"""Cumulative controlled adaptive-basis release campaign for v0.25.3."""

from dataclasses import asdict, dataclass, replace
import hashlib
import json
from types import SimpleNamespace

import numpy as np

from .adaptive_multigaussian_tdvp_v252 import (
    build_adaptive_variational_metric_system_v252,
    quadratic_spin_hamiltonian_from_provider_v252,
)
from .controlled_basis_adaptation_v253 import (
    ControlledBasisSettingsV253,
    SpawnCandidateV253,
    adapt_basis_once_v253,
    coefficient_activation_implicit_step_v253,
    evaluate_spawn_candidate_v253,
    project_adaptive_state_v253,
)
from .controlled_basis_validation_v253 import (
    run_controlled_basis_validation_evidence_v253,
)
from .v252_benchmark import run_v0252_release_benchmark


@dataclass(frozen=True)
class V253AcceptanceThresholds:
    expected_inherited_gates: int = 630
    expected_validation_gates: int = 60
    expected_core_gates: int = 25
    expected_new_gates: int = 85
    expected_total_gates: int = 715


def _raises_v253(callable_object, exceptions, text=None):
    try:
        callable_object()
    except exceptions as exc:
        return text is None or text in str(exc)
    return False


def _core_controls_v253(evidence):
    settings = ControlledBasisSettingsV253().validate()
    event = evidence.odd_spawn_event
    state = event.before
    model = event.model
    activation = evidence.activation_step
    duplicate = SpawnCandidateV253(
        state.q[0], state.p[0], state.widths[0], state.chirps[0], 0, "external"
    )
    static_provider = SimpleNamespace(
        evaluate_snapshot=lambda q: SimpleNamespace(matrices=object())
    )
    canonical = evidence.as_dict()
    payload = json.dumps(
        canonical, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("utf-8")
    digest = hashlib.sha256(payload).hexdigest()
    controls = {
        "multidimensional_request_is_rejected": _raises_v253(
            lambda: replace(settings, multidimensional_nuclear_motion=True).validate(),
            (ValueError,), "does not admit",
        ),
        "full_width_matrix_request_is_rejected": _raises_v253(
            lambda: replace(settings, full_width_matrices=True).validate(),
            (ValueError,), "does not admit",
        ),
        "moving_frame_request_is_rejected": _raises_v253(
            lambda: replace(settings, coordinate_dependent_electronic_frame=True).validate(),
            (ValueError,), "does not admit",
        ),
        "real_molecular_soc_request_is_rejected": _raises_v253(
            lambda: replace(settings, real_molecular_soc_provider=True).validate(),
            (ValueError,), "does not admit",
        ),
        "general_aims_request_is_rejected": _raises_v253(
            lambda: replace(settings, general_aims_branching=True).validate(),
            (ValueError,), "does not admit",
        ),
        "multiple_events_request_is_rejected": _raises_v253(
            lambda: replace(settings, one_event_per_checkpoint=False).validate(),
            (ValueError,), "one topology event",
        ),
        "disabled_lifecycle_procedure_is_rejected": _raises_v253(
            lambda: replace(settings, pruning=False).validate(),
            (ValueError,), "spawn, prune, and merge",
        ),
        "invalid_novelty_gate_is_rejected": _raises_v253(
            lambda: replace(settings, minimum_spawn_novelty=1.0).validate(),
            (ValueError,), "smaller than one",
        ),
        "invalid_activation_age_order_is_rejected": _raises_v253(
            lambda: replace(
                settings,
                minimum_packet_age_steps=4,
                maximum_activation_age_steps=4,
            ).validate(),
            (ValueError,), "must exceed",
        ),
        "static_provider_is_not_lifecycle_admitted": _raises_v253(
            lambda: quadratic_spin_hamiltonian_from_provider_v252(static_provider),
            (TypeError,), "explicit operator provenance",
        ),
        "duplicate_projection_basis_is_rejected": _raises_v253(
            lambda: project_adaptive_state_v253(
                state, model,
                q=[0.0, 0.0], p=[0.0, 0.0], widths=[2.0, 2.0],
                chirps=[0.0, 0.0], event_kind="prune", settings=settings,
            ),
            (ValueError,), "rank deficient",
        ),
        "duplicate_spawn_candidate_is_rejected": (
            evaluate_spawn_candidate_v253(state, model, duplicate, settings=settings).admitted is False
        ),
        "newborn_full_metric_failure_is_exposed": _raises_v253(
            lambda: build_adaptive_variational_metric_system_v252(activation.end, model),
            (ValueError,), "incompatible with its null space",
        ),
        "all_active_activation_call_is_rejected": _raises_v253(
            lambda: coefficient_activation_implicit_step_v253(
                state, model, 0.02, np.asarray([True, True]),
            ),
            (ValueError,), "requires at least one dormant shape",
        ),
        "activation_nonconvergence_is_rejected": _raises_v253(
            lambda: coefficient_activation_implicit_step_v253(
                event.after,
                model,
                0.5,
                np.asarray([True, True, False]),
                settings=replace(
                    settings.tdvp_settings,
                    nonlinear_max_function_evaluations=1,
                ).validate(),
            ),
            (RuntimeError,), "activation implicit midpoint solve failed",
        ),
        "negative_packet_age_is_rejected": _raises_v253(
            lambda: adapt_basis_once_v253(
                state, model, packet_ids=("a", "b"), packet_ages=(2, -1),
                next_packet_serial=2, settings=settings,
            ),
            (ValueError,), "packet age",
        ),
        "duplicate_packet_id_is_rejected": _raises_v253(
            lambda: adapt_basis_once_v253(
                state, model, packet_ids=("a", "a"), packet_ages=(2, 2),
                next_packet_serial=2, settings=settings,
            ),
            (ValueError,), "unique",
        ),
        "packet_serial_collision_is_rejected": _raises_v253(
            lambda: adapt_basis_once_v253(
                state, model, packet_ids=("g000002", "b"), packet_ages=(2, 2),
                next_packet_serial=2, settings=settings,
            ),
            (ValueError,), "collides",
        ),
        "maximum_packet_count_fails_closed": (
            evidence.maximum_count_event.event_kind == "none"
            and evidence.maximum_count_event.reason == "maximum packet count reached"
        ),
        "tampered_projection_loss_is_rejected": _raises_v253(
            lambda: replace(
                event.projection,
                relative_projection_loss=event.projection.relative_projection_loss + 1.0e-3,
            ).validate(),
            (ValueError,), "relative_projection_loss",
        ),
        "tampered_event_id_is_rejected": _raises_v253(
            lambda: replace(event, packet_ids_after=("a", "b", "bad")).validate(),
            (ValueError,), "packet ID",
        ),
        "tampered_candidate_score_is_rejected": _raises_v253(
            lambda: replace(
                event.candidate_evaluations[0],
                residual_capture=event.candidate_evaluations[0].residual_capture + 1.0e-3,
            ).validate(),
            (ValueError,), "residual_capture",
        ),
        "tampered_activation_shape_is_rejected": _raises_v253(
            lambda: replace(
                activation,
                end=replace(
                    activation.end,
                    q=activation.end.q + np.asarray([0.0, 0.0, 1.0e-3]),
                ),
            ).validate(),
            (ValueError,), "inactive newborn shape",
        ),
        "tampered_activation_residual_is_rejected": _raises_v253(
            lambda: replace(
                activation,
                nonlinear_residual=activation.nonlinear_residual + 1.0e-3,
            ).validate(),
            (ValueError,), "nonlinear residual",
        ),
        "canonical_sha256_evidence_fingerprint_is_stable": bool(
            len(evidence.fingerprint()) == 64
            and evidence.fingerprint() == digest
            and digest == hashlib.sha256(payload).hexdigest()
            and all(character in "0123456789abcdef" for character in digest)
        ),
    }
    if len(controls) != 25:
        raise AssertionError(
            f"v0.25.3 must define exactly 25 core gates, found {len(controls)}."
        )
    return {name: bool(value) for name, value in controls.items()}


def run_v0253_release_benchmark(
    thresholds=V253AcceptanceThresholds(), *, memory_probe_policy="proc_self"
):
    inherited = run_v0252_release_benchmark(memory_probe_policy=memory_probe_policy)
    evidence = run_controlled_basis_validation_evidence_v253()
    if len(evidence.audit.checks) != thresholds.expected_validation_gates:
        raise AssertionError("v0.25.3 validation evidence must define 60 gates.")
    core = _core_controls_v253(evidence)
    inherited_checks = {
        f"inherited_v0252::{name}": bool(value)
        for name, value in inherited["acceptance"]["checks"].items()
    }
    validation_checks = {
        f"controlled_basis_validation::{name}": bool(value)
        for name, value in evidence.audit.checks.items()
    }
    core_checks = {
        f"controlled_basis_core::{name}": bool(value)
        for name, value in core.items()
    }
    if len(inherited_checks) != thresholds.expected_inherited_gates:
        raise AssertionError("v0.25.3 must inherit exactly 630 v0.25.2 gates.")
    new_checks = {**validation_checks, **core_checks}
    if len(new_checks) != thresholds.expected_new_gates:
        raise AssertionError("v0.25.3 must define exactly 85 new gates.")
    checks = {**inherited_checks, **new_checks}
    if len(checks) != thresholds.expected_total_gates:
        raise AssertionError("v0.25.3 must define exactly 715 cumulative gates.")
    return {
        "release": "v0.25.3",
        "theme": (
            "controlled residual-driven spawning, coefficient-only newborn activation, "
            "and projection-guarded merge/prune lifecycle"
        ),
        "controlled_basis_validation_evidence": evidence.as_dict(),
        "controlled_basis_core_controls": core,
        "claims": evidence.claims,
        "inherited_v0252": inherited,
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


def save_v0253_release_benchmark(path, **kwargs):
    from pathlib import Path

    result = run_v0253_release_benchmark(**kwargs)
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return target
