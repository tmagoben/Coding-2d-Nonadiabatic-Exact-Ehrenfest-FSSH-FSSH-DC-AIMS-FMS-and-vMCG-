"""Cumulative adaptive-width multi-Gaussian TDVP campaign for v0.25.2."""

from dataclasses import asdict, dataclass, replace
import hashlib
import json
from types import SimpleNamespace

import numpy as np

from .adaptive_multigaussian_tdvp_v252 import (
    AdaptiveVariationalSettingsV252,
    adaptive_implicit_midpoint_tdvp_step_v252,
    build_adaptive_variational_metric_system_v252,
    quadratic_spin_hamiltonian_from_provider_v252,
)
from .adaptive_multigaussian_tdvp_validation_v252 import (
    run_adaptive_multigaussian_validation_evidence_v252,
)
from .multigaussian_tdvp_v251 import solve_variational_metric_v251
from .v251_benchmark import run_v0251_release_benchmark


@dataclass(frozen=True)
class V252AcceptanceThresholds:
    expected_inherited_gates: int = 535
    expected_validation_gates: int = 70
    expected_core_gates: int = 25
    expected_new_gates: int = 95
    expected_total_gates: int = 630


def _raises_v252(callable_object, exceptions, text=None):
    try:
        callable_object()
    except exceptions as exc:
        return text is None or text in str(exc)
    return False


def _core_controls_v252(evidence):
    settings = AdaptiveVariationalSettingsV252().validate()
    receipt = evidence.odd_trajectory.steps[0]
    state = receipt.start
    model = receipt.model
    count = state.parameter_count
    broken_metric = replace(
        receipt.midpoint_system,
        metric=receipt.midpoint_system.metric + 1.0e-3 * np.eye(count),
    )
    broken_rhs = replace(
        receipt.midpoint_system,
        rhs=receipt.midpoint_system.rhs + 1.0e-3,
    )
    broken_velocity = replace(
        receipt.midpoint_system,
        velocity=receipt.midpoint_system.velocity + 1.0e-3,
    )
    broken_solve = replace(
        receipt.midpoint_system.solve_receipt,
        singular_values=0.99
        * receipt.midpoint_system.solve_receipt.singular_values,
    )
    broken_spectrum = replace(receipt.midpoint_system, solve_receipt=broken_solve)
    canonical = evidence.as_dict()
    canonical_payload = json.dumps(
        canonical,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
    ).encode("utf-8")
    canonical_digest = hashlib.sha256(canonical_payload).hexdigest()
    static_provider = SimpleNamespace(
        evaluate_snapshot=lambda q: SimpleNamespace(matrices=object())
    )

    controls = {
        "adaptive_width_disable_request_is_rejected": _raises_v252(
            lambda: replace(settings, adaptive_gaussian_widths=False).validate(),
            (ValueError,),
            "adaptive Gaussian widths as enabled",
        ),
        "spawning_request_is_rejected": _raises_v252(
            lambda: replace(settings, spawning=True).validate(),
            (ValueError,),
            "spawning",
        ),
        "pruning_request_is_rejected": _raises_v252(
            lambda: replace(settings, pruning=True).validate(),
            (ValueError,),
            "pruning",
        ),
        "coordinate_dependent_frame_request_is_rejected": _raises_v252(
            lambda: replace(
                settings, coordinate_dependent_electronic_frame=True
            ).validate(),
            (ValueError,),
            "coordinate_dependent_electronic_frame",
        ),
        "multidimensional_request_is_rejected": _raises_v252(
            lambda: replace(
                settings, multidimensional_nuclear_motion=True
            ).validate(),
            (ValueError,),
            "multidimensional_nuclear_motion",
        ),
        "full_width_matrix_request_is_rejected": _raises_v252(
            lambda: replace(settings, full_width_matrices=True).validate(),
            (ValueError,),
            "full_width_matrices",
        ),
        "real_molecular_provider_request_is_rejected": _raises_v252(
            lambda: replace(settings, real_molecular_soc_provider=True).validate(),
            (ValueError,),
            "real_molecular_soc_provider",
        ),
        "non_svd_metric_solver_is_rejected": _raises_v252(
            lambda: replace(settings, metric_solver="normal equations").validate(),
            (ValueError,),
            "SVD metric solver is frozen",
        ),
        "non_midpoint_integrator_is_rejected": _raises_v252(
            lambda: replace(settings, integrator="velocity Verlet").validate(),
            (ValueError,),
            "implicit integrator is frozen",
        ),
        "wrong_width_coordinates_are_rejected": _raises_v252(
            lambda: replace(settings, width_coordinates="raw alpha only").validate(),
            (ValueError,),
            "width coordinates are frozen",
        ),
        "inverted_width_domain_is_rejected": _raises_v252(
            lambda: replace(
                settings, minimum_width=2.0, maximum_width=1.0
            ).validate(),
            (ValueError,),
            "maximum_width must exceed minimum_width",
        ),
        "incompatible_null_rhs_is_rejected": _raises_v252(
            lambda: solve_variational_metric_v251(
                np.diag([1.0, 0.0]),
                np.asarray([0.0, 1.0]),
                settings=settings,
            ),
            (ValueError,),
            "incompatible with its null space",
        ),
        "indefinite_metric_is_rejected": _raises_v252(
            lambda: solve_variational_metric_v251(
                np.diag([1.0, -0.1]), np.zeros(2), settings=settings
            ),
            (ValueError,),
            "not positive semidefinite",
        ),
        "static_provider_is_not_adaptive_trajectory_admitted": _raises_v252(
            lambda: quadratic_spin_hamiltonian_from_provider_v252(static_provider),
            (TypeError,),
            "explicit operator provenance",
        ),
        "nonlinear_nonconvergence_is_rejected": _raises_v252(
            lambda: adaptive_implicit_midpoint_tdvp_step_v252(
                state,
                model,
                0.8,
                settings=replace(
                    settings, nonlinear_max_function_evaluations=1
                ).validate(),
            ),
            (RuntimeError,),
            "adaptive implicit midpoint TDVP solve failed",
        ),
        "width_below_configured_domain_is_rejected": _raises_v252(
            lambda: build_adaptive_variational_metric_system_v252(
                replace(state, widths=np.asarray([1.0e-10, state.widths[1]])),
                model,
                settings=settings,
            ),
            (ValueError,),
            "configured minimum",
        ),
        "chirp_above_configured_domain_is_rejected": _raises_v252(
            lambda: build_adaptive_variational_metric_system_v252(
                replace(
                    state,
                    chirps=np.asarray(
                        [settings.maximum_absolute_chirp * 2.0, state.chirps[1]]
                    ),
                ),
                model,
                settings=settings,
            ),
            (ValueError,),
            "configured maximum",
        ),
        "excessive_per_step_log_width_change_is_rejected": _raises_v252(
            lambda: replace(
                receipt,
                settings=replace(
                    settings, maximum_step_log_width_change=1.0e-8
                ).validate(),
            ).validate(),
            (ValueError,),
            "width changed too much",
        ),
        "tampered_endpoint_is_rejected": _raises_v252(
            lambda: replace(
                receipt,
                end=replace(receipt.end, p=receipt.end.p + 1.0e-3),
            ).validate(),
            (ValueError,),
        ),
        "tampered_midpoint_metric_is_rejected": _raises_v252(
            lambda: replace(receipt, midpoint_system=broken_metric).validate(),
            (ValueError,),
        ),
        "tampered_midpoint_rhs_is_rejected": _raises_v252(
            lambda: replace(receipt, midpoint_system=broken_rhs).validate(),
            (ValueError,),
        ),
        "tampered_midpoint_velocity_is_rejected": _raises_v252(
            lambda: replace(receipt, midpoint_system=broken_velocity).validate(),
            (ValueError,),
        ),
        "tampered_metric_spectrum_is_rejected": _raises_v252(
            lambda: replace(receipt, midpoint_system=broken_spectrum).validate(),
            (ValueError,),
        ),
        "tampered_nonlinear_residual_is_rejected": _raises_v252(
            lambda: replace(
                receipt,
                nonlinear_residual=receipt.nonlinear_residual + 1.0e-3,
            ).validate(),
            (ValueError,),
            "adaptive nonlinear residual",
        ),
        "canonical_sha256_evidence_fingerprint_is_stable": bool(
            len(evidence.fingerprint()) == 64
            and evidence.fingerprint() == canonical_digest
            and canonical_digest
            == hashlib.sha256(canonical_payload).hexdigest()
            and all(
                character in "0123456789abcdef"
                for character in evidence.fingerprint()
            )
        ),
    }
    if len(controls) != 25:
        raise AssertionError(
            f"v0.25.2 must define exactly 25 core gates, found {len(controls)}."
        )
    return {name: bool(value) for name, value in controls.items()}


def run_v0252_release_benchmark(
    thresholds=V252AcceptanceThresholds(),
    *,
    memory_probe_policy="proc_self",
):
    inherited = run_v0251_release_benchmark(
        memory_probe_policy=memory_probe_policy
    )
    evidence = run_adaptive_multigaussian_validation_evidence_v252()
    if len(evidence.audit.checks) != thresholds.expected_validation_gates:
        raise AssertionError("v0.25.2 validation evidence must define 70 gates.")
    core = _core_controls_v252(evidence)
    if len(core) != thresholds.expected_core_gates:
        raise AssertionError("v0.25.2 must define exactly 25 core gates.")

    inherited_checks = {
        f"inherited_v0251::{name}": bool(value)
        for name, value in inherited["acceptance"]["checks"].items()
    }
    validation_checks = {
        f"adaptive_multigaussian_validation::{name}": bool(value)
        for name, value in evidence.audit.checks.items()
    }
    core_checks = {
        f"adaptive_multigaussian_core::{name}": bool(value)
        for name, value in core.items()
    }
    if len(inherited_checks) != thresholds.expected_inherited_gates:
        raise AssertionError("v0.25.2 must inherit exactly 535 v0.25.1 gates.")
    new_checks = {**validation_checks, **core_checks}
    if len(new_checks) != thresholds.expected_new_gates:
        raise AssertionError("v0.25.2 must define exactly 95 new gates.")
    checks = {**inherited_checks, **new_checks}
    if len(checks) != thresholds.expected_total_gates:
        raise AssertionError("v0.25.2 must define exactly 630 cumulative gates.")

    return {
        "release": "v0.25.2",
        "theme": (
            "one-dimensional adaptive log-width/quadratic-chirp multi-Gaussian "
            "McLachlan TDVP with compatible-null full SVD and implicit midpoint"
        ),
        "adaptive_multigaussian_validation_evidence": evidence.as_dict(),
        "adaptive_multigaussian_core_controls": core,
        "claims": evidence.claims,
        "inherited_v0251": inherited,
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

