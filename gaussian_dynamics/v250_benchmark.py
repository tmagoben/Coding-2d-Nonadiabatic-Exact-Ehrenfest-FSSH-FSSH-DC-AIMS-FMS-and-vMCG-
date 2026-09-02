"""Cumulative symmetric variational SOC acceptance campaign for v0.25.0."""

from dataclasses import asdict, dataclass, replace
import hashlib
import json
from types import SimpleNamespace

import numpy as np

from .variational_soc_dynamics_v250 import (
    CanonicalVariationalSOCStateV250,
    VariationalSOCIntegratorSettingsV250,
    symmetric_variational_soc_step_v250,
)
from .variational_soc_validation_v250 import (
    _ScaledOverlapProviderV250,
    run_variational_soc_validation_evidence_v250,
)
from .v242_benchmark import run_v0242_release_benchmark


@dataclass(frozen=True)
class V250AcceptanceThresholds:
    expected_inherited_gates: int = 400
    expected_validation_gates: int = 45
    expected_core_gates: int = 15
    expected_new_gates: int = 60
    expected_total_gates: int = 460


def _raises_v250(callable_object, exceptions, text=None):
    try:
        callable_object()
    except exceptions as exc:
        return text is None or text in str(exc)
    return False


def _core_controls_v250(evidence):
    settings = VariationalSOCIntegratorSettingsV250().validate()
    receipt = evidence.odd_trajectory.steps[0]
    initial = receipt.start

    expansion_rejected = _raises_v250(
        lambda: symmetric_variational_soc_step_v250(
            initial, _ScaledOverlapProviderV250(1.01), 0.2
        ),
        (ValueError,),
        "physically inconsistent",
    )
    rank_loss_rejected = _raises_v250(
        lambda: symmetric_variational_soc_step_v250(
            initial, _ScaledOverlapProviderV250(0.5), 0.2
        ),
        (ValueError,),
        "not trajectory ready",
    )
    static_provider = SimpleNamespace(
        evaluate_snapshot=lambda q: SimpleNamespace(matrices=object()),
        snapshot_overlap=lambda left, right: np.eye(4),
    )

    phase_rotated_end = replace(
        receipt.end,
        electronic_coefficients=(
            np.exp(1.0e-3j) * receipt.end.electronic_coefficients
        ),
    )
    broken_metrics = dict(receipt.transport_metrics)
    broken_metrics["minimum_singular_value"] = 0.2
    canonical_payload = json.dumps(
        evidence.as_dict(),
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
    ).encode("utf-8")
    canonical_digest = hashlib.sha256(canonical_payload).hexdigest()

    return {
        "full_multi_gaussian_tdvp_request_is_rejected": _raises_v250(
            lambda: replace(settings, full_multi_gaussian_tdvp=True).validate(),
            (ValueError,),
            "does not admit full multi-Gaussian TDVP",
        ),
        "adaptive_width_request_is_rejected": _raises_v250(
            lambda: replace(settings, adaptive_gaussian_widths=True).validate(),
            (ValueError,),
            "adaptive Gaussian-width",
        ),
        "coordinate_dependent_mass_verlet_request_is_rejected": _raises_v250(
            lambda: replace(settings, coordinate_dependent_mass=True).validate(),
            (ValueError,),
            "coordinate-dependent generalized mass",
        ),
        "non_svd_polar_algorithm_is_rejected": _raises_v250(
            lambda: replace(settings, polar_algorithm="raw overlap").validate(),
            (ValueError,),
            "SVD-polar algorithm is frozen",
        ),
        "spectrally_expansive_overlap_is_rejected": expansion_rejected,
        "rank_lost_overlap_is_rejected": rank_loss_rejected,
        "static_soc_snapshot_is_not_trajectory_admitted": _raises_v250(
            lambda: symmetric_variational_soc_step_v250(
                initial, static_provider, 0.2
            ),
            (TypeError,),
            "full H, K, D, and mass",
        ),
        "tampered_endpoint_momentum_is_rejected": _raises_v250(
            lambda: replace(
                receipt,
                end=replace(receipt.end, p=receipt.end.p + 1.0e-3),
            ).validate(),
            (ValueError,),
            "endpoint momentum disagrees",
        ),
        "tampered_endpoint_spinor_is_rejected": _raises_v250(
            lambda: replace(receipt, end=phase_rotated_end).validate(),
            (ValueError,),
            "electronic endpoint disagrees",
        ),
        "tampered_polar_transport_is_rejected": _raises_v250(
            lambda: replace(
                receipt,
                transport_end_to_start=np.diag([1.0, 1.0, 1.0, -1.0]),
            ).validate(),
            (ValueError,),
            "polar transport disagrees",
        ),
        "tampered_singular_values_are_rejected": _raises_v250(
            lambda: replace(
                receipt, singular_values=0.95 * receipt.singular_values
            ).validate(),
            (ValueError,),
            "singular values disagree",
        ),
        "tampered_transport_metrics_are_rejected": _raises_v250(
            lambda: replace(receipt, transport_metrics=broken_metrics).validate(),
            (ValueError,),
            "metric minimum_singular_value",
        ),
        "tampered_endpoint_mass_is_rejected": _raises_v250(
            lambda: replace(
                receipt,
                mass_matrix_end_au=receipt.mass_matrix_end_au + 1.0e-3,
            ).validate(),
            (ValueError,),
            "coordinate-dependent mass",
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
            and canonical_digest == hashlib.sha256(canonical_payload).hexdigest()
        ),
    }


def run_v0250_release_benchmark(
    thresholds=V250AcceptanceThresholds(),
    *,
    memory_probe_policy="proc_self",
):
    inherited = run_v0242_release_benchmark(
        memory_probe_policy=memory_probe_policy
    )
    evidence = run_variational_soc_validation_evidence_v250()
    if len(evidence.audit.checks) != thresholds.expected_validation_gates:
        raise AssertionError("v0.25.0 validation evidence must define exactly 45 gates.")
    core = _core_controls_v250(evidence)
    if len(core) != thresholds.expected_core_gates:
        raise AssertionError("v0.25.0 must define exactly 15 core gates.")

    inherited_checks = {
        f"inherited_v0242::{name}": bool(value)
        for name, value in inherited["acceptance"]["checks"].items()
    }
    validation_checks = {
        f"variational_soc_validation::{name}": bool(value)
        for name, value in evidence.audit.checks.items()
    }
    core_checks = {
        f"variational_soc_core::{name}": bool(value)
        for name, value in core.items()
    }
    if len(inherited_checks) != thresholds.expected_inherited_gates:
        raise AssertionError("v0.25.0 must inherit exactly 400 v0.24.2 gates.")
    new_checks = {**validation_checks, **core_checks}
    if len(new_checks) != thresholds.expected_new_gates:
        raise AssertionError("v0.25.0 must define exactly 60 new gates.")
    checks = {**inherited_checks, **new_checks}
    if len(checks) != thresholds.expected_total_gates:
        raise AssertionError("v0.25.0 must define exactly 460 cumulative gates.")

    return {
        "release": "v0.25.0",
        "theme": (
            "restricted time-dependent-variational SOC dynamics with symmetric "
            "constant-mass Verlet, endpoint Strang propagation, and SVD-polar transport"
        ),
        "variational_soc_validation_evidence": evidence.as_dict(),
        "variational_soc_core_controls": core,
        "claims": evidence.claims,
        "inherited_v0242": inherited,
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
