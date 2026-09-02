"""Cumulative frozen-width multi-Gaussian TDVP acceptance campaign for v0.25.1."""

from dataclasses import asdict, dataclass, replace
import hashlib
import json
from types import SimpleNamespace

import numpy as np

from .multigaussian_tdvp_v251 import (
    VariationalMetricSettingsV251,
    implicit_midpoint_tdvp_step_v251,
    quadratic_spin_hamiltonian_from_provider_v251,
    solve_variational_metric_v251,
)
from .multigaussian_tdvp_validation_v251 import (
    run_multigaussian_tdvp_validation_evidence_v251,
)
from .v250_benchmark import run_v0250_release_benchmark


@dataclass(frozen=True)
class V251AcceptanceThresholds:
    expected_inherited_gates: int = 460
    expected_validation_gates: int = 55
    expected_core_gates: int = 20
    expected_new_gates: int = 75
    expected_total_gates: int = 535


def _raises_v251(callable_object, exceptions, text=None):
    try:
        callable_object()
    except exceptions as exc:
        return text is None or text in str(exc)
    return False


def _core_controls_v251(evidence):
    settings = VariationalMetricSettingsV251().validate()
    receipt = evidence.odd_trajectory.steps[0]
    state = receipt.start
    model = receipt.model
    broken_metric = replace(
        receipt.midpoint_system,
        metric=receipt.midpoint_system.metric
        + 1.0e-3 * np.eye(state.parameter_count),
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
        "adaptive_width_request_is_rejected": _raises_v251(
            lambda: replace(settings, adaptive_gaussian_widths=True).validate(),
            (ValueError,),
            "adaptive_gaussian_widths",
        ),
        "spawning_request_is_rejected": _raises_v251(
            lambda: replace(settings, spawning=True).validate(),
            (ValueError,),
            "spawning",
        ),
        "pruning_request_is_rejected": _raises_v251(
            lambda: replace(settings, pruning=True).validate(),
            (ValueError,),
            "pruning",
        ),
        "coordinate_dependent_frame_request_is_rejected": _raises_v251(
            lambda: replace(
                settings, coordinate_dependent_electronic_frame=True
            ).validate(),
            (ValueError,),
            "coordinate_dependent_electronic_frame",
        ),
        "multidimensional_request_is_rejected": _raises_v251(
            lambda: replace(
                settings, multidimensional_nuclear_motion=True
            ).validate(),
            (ValueError,),
            "multidimensional_nuclear_motion",
        ),
        "real_molecular_provider_request_is_rejected": _raises_v251(
            lambda: replace(settings, real_molecular_soc_provider=True).validate(),
            (ValueError,),
            "real_molecular_soc_provider",
        ),
        "non_svd_metric_solver_is_rejected": _raises_v251(
            lambda: replace(settings, metric_solver="normal equations").validate(),
            (ValueError,),
            "SVD metric solver is frozen",
        ),
        "non_midpoint_integrator_is_rejected": _raises_v251(
            lambda: replace(settings, integrator="velocity Verlet").validate(),
            (ValueError,),
            "implicit integrator is frozen",
        ),
        "incompatible_null_rhs_is_rejected": _raises_v251(
            lambda: solve_variational_metric_v251(
                np.diag([1.0, 0.0]), np.asarray([0.0, 1.0])
            ),
            (ValueError,),
            "incompatible with its null space",
        ),
        "indefinite_metric_is_rejected": _raises_v251(
            lambda: solve_variational_metric_v251(
                np.diag([1.0, -0.1]), np.zeros(2)
            ),
            (ValueError,),
            "not positive semidefinite",
        ),
        "static_provider_is_not_quadratic_trajectory_admitted": _raises_v251(
            lambda: quadratic_spin_hamiltonian_from_provider_v251(static_provider),
            (TypeError,),
            "explicit operator provenance",
        ),
        "nonlinear_nonconvergence_is_rejected": _raises_v251(
            lambda: implicit_midpoint_tdvp_step_v251(
                state,
                model,
                0.8,
                settings=replace(
                    settings, nonlinear_max_function_evaluations=1
                ).validate(),
            ),
            (RuntimeError,),
            "implicit midpoint TDVP solve failed",
        ),
        "tampered_endpoint_is_rejected": _raises_v251(
            lambda: replace(
                receipt,
                end=replace(receipt.end, p=receipt.end.p + 1.0e-3),
            ).validate(),
            (ValueError,),
        ),
        "tampered_midpoint_metric_is_rejected": _raises_v251(
            lambda: replace(receipt, midpoint_system=broken_metric).validate(),
            (ValueError,),
        ),
        "tampered_midpoint_rhs_is_rejected": _raises_v251(
            lambda: replace(receipt, midpoint_system=broken_rhs).validate(),
            (ValueError,),
        ),
        "tampered_midpoint_velocity_is_rejected": _raises_v251(
            lambda: replace(receipt, midpoint_system=broken_velocity).validate(),
            (ValueError,),
        ),
        "tampered_metric_spectrum_is_rejected": _raises_v251(
            lambda: replace(receipt, midpoint_system=broken_spectrum).validate(),
            (ValueError,),
        ),
        "tampered_nonlinear_residual_is_rejected": _raises_v251(
            lambda: replace(
                receipt,
                nonlinear_residual=receipt.nonlinear_residual + 1.0e-3,
            ).validate(),
            (ValueError,),
            "nonlinear midpoint residual",
        ),
        "validation_evidence_fingerprint_is_sha256": bool(
            len(evidence.fingerprint()) == 64
            and all(
                character in "0123456789abcdef"
                for character in evidence.fingerprint()
            )
        ),
        "canonical_evidence_serialization_is_stable": bool(
            evidence.fingerprint() == canonical_digest
            and canonical_digest
            == hashlib.sha256(canonical_payload).hexdigest()
        ),
    }
    if len(controls) != 20:
        raise AssertionError("v0.25.1 must define exactly 20 core gates.")
    return {name: bool(value) for name, value in controls.items()}


def run_v0251_release_benchmark(
    thresholds=V251AcceptanceThresholds(),
    *,
    memory_probe_policy="proc_self",
):
    inherited = run_v0250_release_benchmark(
        memory_probe_policy=memory_probe_policy
    )
    evidence = run_multigaussian_tdvp_validation_evidence_v251()
    if len(evidence.audit.checks) != thresholds.expected_validation_gates:
        raise AssertionError("v0.25.1 validation evidence must define 55 gates.")
    core = _core_controls_v251(evidence)
    if len(core) != thresholds.expected_core_gates:
        raise AssertionError("v0.25.1 must define exactly 20 core gates.")

    inherited_checks = {
        f"inherited_v0250::{name}": bool(value)
        for name, value in inherited["acceptance"]["checks"].items()
    }
    validation_checks = {
        f"multigaussian_tdvp_validation::{name}": bool(value)
        for name, value in evidence.audit.checks.items()
    }
    core_checks = {
        f"multigaussian_tdvp_core::{name}": bool(value)
        for name, value in core.items()
    }
    if len(inherited_checks) != thresholds.expected_inherited_gates:
        raise AssertionError("v0.25.1 must inherit exactly 460 v0.25.0 gates.")
    new_checks = {**validation_checks, **core_checks}
    if len(new_checks) != thresholds.expected_new_gates:
        raise AssertionError("v0.25.1 must define exactly 75 new gates.")
    checks = {**inherited_checks, **new_checks}
    if len(checks) != thresholds.expected_total_gates:
        raise AssertionError("v0.25.1 must define exactly 535 cumulative gates.")

    return {
        "release": "v0.25.1",
        "theme": (
            "one-dimensional frozen-width multi-Gaussian McLachlan TDVP with "
            "a compatible-null SVD metric solve and fully implicit midpoint"
        ),
        "multigaussian_tdvp_validation_evidence": evidence.as_dict(),
        "multigaussian_tdvp_core_controls": core,
        "claims": evidence.claims,
        "inherited_v0250": inherited,
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
