"""Deterministic validation evidence for the v0.25.0 variational SOC integrator."""

from dataclasses import dataclass

import numpy as np

from .analytic_soc_models_v220 import (
    AnalyticDoubletSOCProviderV220,
    AnalyticSingletTripletSOCProviderV220,
    SingletTripletSOCConfigV220,
)
from .complex_gauge_v21 import (
    GaugeTransformedOperatorProviderV21,
    PhaseMixingGaugeV21,
    random_unitary_v21,
)
from .variational_soc_dynamics_v250 import (
    ELECTRONIC_INTEGRATOR_V250,
    GENERAL_TDVP_INTEGRATOR_V250,
    POLAR_ALGORITHM_V250,
    RESTRICTED_NUCLEAR_INTEGRATOR_V250,
    RESTRICTED_TDVP_ANSATZ_V250,
    V250_TRAJECTORY_CLAIMS,
    VARIATIONAL_SOC_SCHEMA_V250,
    CanonicalVariationalSOCStateV250,
    _sha256_v250,
    reverse_variational_soc_trajectory_v250,
    run_symmetric_variational_soc_dynamics_v250,
)


VARIATIONAL_SOC_VALIDATION_SCHEMA_V250 = (
    "gnd-symmetric-variational-soc-validation-v0.25.0"
)
V250_CONVERGENCE_DT_AU = (0.8, 0.4, 0.2, 0.1)
V250_CONVERGENCE_FINAL_TIME_AU = 20.0


def _initial_state_v250():
    coefficients = np.asarray(
        [0.67 + 0.11j, 0.19 - 0.28j, 0.41 + 0.17j, -0.09j],
        dtype=complex,
    )
    coefficients /= np.linalg.norm(coefficients)
    return CanonicalVariationalSOCStateV250(
        np.asarray([-0.7]), np.asarray([8.0]), coefficients
    ).validate()


def _gauge_v250():
    return PhaseMixingGaugeV21(
        random_unitary_v21(4, 25001),
        np.asarray([[0.17], [-0.11], [0.23], [-0.07]]),
        np.asarray([0.20, -0.30, 0.10, 0.40]),
    )


class _ScaledOverlapProviderV250:
    def __init__(self, singular_value=0.97):
        self.base = AnalyticDoubletSOCProviderV220()
        self.singular_value = float(singular_value)

    def evaluate_snapshot(self, q):
        return self.base.evaluate_snapshot(q)

    def snapshot_overlap(self, left, right):
        return self.singular_value * np.eye(4, dtype=complex)


def _phase_aligned_state_distance_v250(left, right):
    overlap = np.vdot(
        right.electronic_coefficients, left.electronic_coefficients
    )
    phase = 1.0 + 0.0j if abs(overlap) < 1.0e-30 else overlap / abs(overlap)
    electronic = left.electronic_coefficients / phase
    return float(
        np.sqrt(
            np.linalg.norm(left.q - right.q) ** 2
            + np.linalg.norm(left.p - right.p) ** 2
            + np.linalg.norm(electronic - right.electronic_coefficients) ** 2
        )
    )


def _trajectory_finite_v250(trajectory):
    states = [trajectory.initial_state] + [step.end for step in trajectory.steps]
    return bool(
        all(
            np.all(np.isfinite(state.q))
            and np.all(np.isfinite(state.p))
            and np.all(np.isfinite(state.electronic_coefficients))
            for state in states
        )
    )


@dataclass(frozen=True)
class VariationalSOCValidationAuditV250:
    checks: dict
    metrics: dict
    thresholds: dict
    passed: bool

    def validate(self):
        if not isinstance(self.checks, dict) or len(self.checks) != 45:
            raise ValueError("v0.25.0 variational SOC audit requires exactly 45 gates.")
        if any(type(value) is not bool for value in self.checks.values()):
            raise TypeError("every v0.25.0 variational SOC gate must be Boolean.")
        if type(self.passed) is not bool or self.passed != bool(all(self.checks.values())):
            raise ValueError("v0.25.0 variational SOC audit result is inconsistent.")
        _sha256_v250(self.metrics)
        _sha256_v250(self.thresholds)
        return self

    def as_dict(self):
        self.validate()
        return {
            "checks": dict(self.checks),
            "metrics": self.metrics,
            "thresholds": self.thresholds,
            "passed": bool(self.passed),
        }


@dataclass(frozen=True)
class VariationalSOCValidationEvidenceV250:
    even_trajectory: object
    even_reverse: object
    odd_trajectory: object
    odd_reverse: object
    gauge_trajectory: object
    contraction_trajectory: object
    zero_soc_enabled_trajectory: object
    zero_soc_disabled_trajectory: object
    convergence_trajectories: tuple
    audit: VariationalSOCValidationAuditV250

    def validate(self):
        for trajectory in (
            self.even_trajectory,
            self.even_reverse,
            self.odd_trajectory,
            self.odd_reverse,
            self.gauge_trajectory,
            self.contraction_trajectory,
            self.zero_soc_enabled_trajectory,
            self.zero_soc_disabled_trajectory,
            *self.convergence_trajectories,
        ):
            trajectory.validate()
        self.audit.validate()
        if not self.audit.passed:
            failed = ", ".join(
                name for name, passed in self.audit.checks.items() if not passed
            )
            raise ValueError("v0.25.0 variational SOC evidence failed: " + failed)
        return self

    @property
    def claims(self):
        return dict(V250_TRAJECTORY_CLAIMS)

    @staticmethod
    def _convergence_receipt(trajectory):
        initial_energy = (
            trajectory.steps[0].energy_start_hartree
            if trajectory.steps
            else 0.0
        )
        return {
            "dt_au": (
                float(trajectory.steps[0].dt_au) if trajectory.steps else 0.0
            ),
            "step_count": len(trajectory.steps),
            "final_state": trajectory.final_state.as_dict(),
            "initial_energy_hartree": float(initial_energy),
            "endpoint_energies_hartree": [
                float(step.energy_end_hartree) for step in trajectory.steps
            ],
            "maximum_norm_drift": trajectory.maximum_norm_drift,
            "fingerprint": trajectory.fingerprint(),
        }

    def as_dict(self):
        self.validate()
        return {
            "schema": VARIATIONAL_SOC_VALIDATION_SCHEMA_V250,
            "trajectory_schema": VARIATIONAL_SOC_SCHEMA_V250,
            "decisions": {
                "validated_ansatz": RESTRICTED_TDVP_ANSATZ_V250,
                "restricted_nuclear_integrator": RESTRICTED_NUCLEAR_INTEGRATOR_V250,
                "general_tdvp_integrator": GENERAL_TDVP_INTEGRATOR_V250,
                "electronic_integrator": ELECTRONIC_INTEGRATOR_V250,
                "polar_algorithm": POLAR_ALGORITHM_V250,
            },
            "canonical_trajectories": {
                "even": self.even_trajectory.as_dict(),
                "even_reverse": self.even_reverse.as_dict(),
                "odd": self.odd_trajectory.as_dict(),
                "odd_reverse": self.odd_reverse.as_dict(),
                "odd_complex_gauge": self.gauge_trajectory.as_dict(),
                "contractive_overlap": self.contraction_trajectory.as_dict(),
                "zero_soc_enabled": self.zero_soc_enabled_trajectory.as_dict(),
                "zero_soc_disabled": self.zero_soc_disabled_trajectory.as_dict(),
            },
            "convergence_receipts": [
                self._convergence_receipt(trajectory)
                for trajectory in self.convergence_trajectories
            ],
            "audit": self.audit.as_dict(),
            "claims": self.claims,
        }

    def fingerprint(self):
        return _sha256_v250(self.as_dict())


def run_variational_soc_validation_evidence_v250():
    """Run even/odd, gauge, reversibility, contraction, and convergence evidence."""

    initial = _initial_state_v250()
    even_provider = AnalyticSingletTripletSOCProviderV220()
    odd_provider = AnalyticDoubletSOCProviderV220()
    even = run_symmetric_variational_soc_dynamics_v250(
        initial, even_provider, dt_au=0.4, steps=50
    )
    odd = run_symmetric_variational_soc_dynamics_v250(
        initial, odd_provider, dt_au=0.4, steps=50
    )
    even_reverse = reverse_variational_soc_trajectory_v250(
        even, AnalyticSingletTripletSOCProviderV220()
    )
    odd_reverse = reverse_variational_soc_trajectory_v250(
        odd, AnalyticDoubletSOCProviderV220()
    )

    gauge = _gauge_v250()
    transformed_initial = CanonicalVariationalSOCStateV250(
        initial.q,
        initial.p,
        gauge.matrix(initial.q).conj().T @ initial.electronic_coefficients,
    ).validate()
    gauge_trajectory = run_symmetric_variational_soc_dynamics_v250(
        transformed_initial,
        GaugeTransformedOperatorProviderV21(
            AnalyticDoubletSOCProviderV220(), gauge
        ),
        dt_au=0.4,
        steps=50,
    )

    contraction = run_symmetric_variational_soc_dynamics_v250(
        initial, _ScaledOverlapProviderV250(), dt_au=0.4, steps=10
    )
    zero_enabled = run_symmetric_variational_soc_dynamics_v250(
        initial,
        AnalyticSingletTripletSOCProviderV220(
            SingletTripletSOCConfigV220(soc_scale=0.0, soc_enabled=True)
        ),
        dt_au=0.4,
        steps=20,
    )
    zero_disabled = run_symmetric_variational_soc_dynamics_v250(
        initial,
        AnalyticSingletTripletSOCProviderV220(
            SingletTripletSOCConfigV220(soc_scale=0.0, soc_enabled=False)
        ),
        dt_au=0.4,
        steps=20,
    )

    convergence = tuple(
        run_symmetric_variational_soc_dynamics_v250(
            initial,
            AnalyticSingletTripletSOCProviderV220(),
            dt_au=dt,
            steps=int(round(V250_CONVERGENCE_FINAL_TIME_AU / dt)),
        )
        for dt in V250_CONVERGENCE_DT_AU
    )
    state_changes = [
        _phase_aligned_state_distance_v250(left.final_state, right.final_state)
        for left, right in zip(convergence[:-1], convergence[1:])
    ]
    state_ratios = [
        state_changes[index + 1] / state_changes[index]
        for index in range(len(state_changes) - 1)
    ]
    energy_drifts = [
        trajectory.maximum_absolute_energy_drift_hartree
        for trajectory in convergence
    ]
    energy_ratios = [
        energy_drifts[index + 1] / energy_drifts[index]
        for index in range(len(energy_drifts) - 1)
    ]

    expected_gauge_coefficients = (
        gauge.matrix(gauge_trajectory.final_state.q).conj().T
        @ odd.final_state.electronic_coefficients
    )
    gauge_coefficient_error = float(
        np.linalg.norm(
            gauge_trajectory.final_state.electronic_coefficients
            - expected_gauge_coefficients
        )
    )
    thresholds = {
        "maximum_norm_drift": 5.0e-13,
        "maximum_energy_drift_hartree": 1.0e-9,
        "maximum_reversibility_error": 5.0e-13,
        "maximum_gauge_error": 5.0e-13,
        "minimum_second_order_ratio": 0.24,
        "maximum_second_order_ratio": 0.26,
        "minimum_contractive_singular_value": 0.9,
        "matrix_residual": 5.0e-12,
    }
    all_canonical_steps = (
        *even.steps,
        *odd.steps,
        *gauge_trajectory.steps,
        *contraction.steps,
    )
    transport_unitarity = [
        float(
            np.linalg.norm(
                step.transport_end_to_start.conj().T
                @ step.transport_end_to_start
                - np.eye(len(step.singular_values)),
                ord="fro",
            )
        )
        for step in all_canonical_steps
    ]
    checks = {
        "validation_schema_is_v0250": VARIATIONAL_SOC_VALIDATION_SCHEMA_V250.endswith("v0.25.0"),
        "restricted_tdvp_ansatz_is_explicit": "single canonical nuclear packet" in RESTRICTED_TDVP_ANSATZ_V250,
        "nuclear_verlet_scope_is_constant_mass_canonical": "constant-mass canonical" in RESTRICTED_NUCLEAR_INTEGRATOR_V250,
        "general_tdvp_recommends_implicit_midpoint": "implicit midpoint" in GENERAL_TDVP_INTEGRATOR_V250,
        "electronic_step_is_endpoint_strang": "Strang" in ELECTRONIC_INTEGRATOR_V250,
        "polar_transport_is_svd_computed": POLAR_ALGORITHM_V250.startswith("SVD"),
        "even_trajectory_is_finite": _trajectory_finite_v250(even),
        "even_electronic_norm_is_preserved": even.maximum_norm_drift <= thresholds["maximum_norm_drift"],
        "even_energy_drift_is_bounded": even.maximum_absolute_energy_drift_hartree <= thresholds["maximum_energy_drift_hartree"],
        "even_reverse_coordinate_returns": np.max(np.abs(even_reverse.final_state.q - initial.q)) <= thresholds["maximum_reversibility_error"],
        "even_reverse_momentum_returns": np.max(np.abs(even_reverse.final_state.p - initial.p)) <= thresholds["maximum_reversibility_error"],
        "even_reverse_spinor_returns": np.max(np.abs(even_reverse.final_state.electronic_coefficients - initial.electronic_coefficients)) <= thresholds["maximum_reversibility_error"],
        "even_reverse_time_returns": abs(even_reverse.final_state.time_au - initial.time_au) <= thresholds["maximum_reversibility_error"],
        "odd_trajectory_is_finite": _trajectory_finite_v250(odd),
        "odd_electronic_norm_is_preserved": odd.maximum_norm_drift <= thresholds["maximum_norm_drift"],
        "odd_energy_drift_is_bounded": odd.maximum_absolute_energy_drift_hartree <= thresholds["maximum_energy_drift_hartree"],
        "odd_reverse_coordinate_returns": np.max(np.abs(odd_reverse.final_state.q - initial.q)) <= thresholds["maximum_reversibility_error"],
        "odd_reverse_momentum_returns": np.max(np.abs(odd_reverse.final_state.p - initial.p)) <= thresholds["maximum_reversibility_error"],
        "odd_reverse_spinor_returns": np.max(np.abs(odd_reverse.final_state.electronic_coefficients - initial.electronic_coefficients)) <= thresholds["maximum_reversibility_error"],
        "odd_reverse_time_returns": abs(odd_reverse.final_state.time_au - initial.time_au) <= thresholds["maximum_reversibility_error"],
        "complex_gauge_coordinate_covariance": np.max(np.abs(gauge_trajectory.final_state.q - odd.final_state.q)) <= thresholds["maximum_gauge_error"],
        "complex_gauge_momentum_covariance": np.max(np.abs(gauge_trajectory.final_state.p - odd.final_state.p)) <= thresholds["maximum_gauge_error"],
        "complex_gauge_spinor_covariance": gauge_coefficient_error <= thresholds["maximum_gauge_error"],
        "complex_gauge_norm_preservation": gauge_trajectory.maximum_norm_drift <= thresholds["maximum_norm_drift"],
        "complex_gauge_energy_drift_matches": abs(gauge_trajectory.maximum_absolute_energy_drift_hartree - odd.maximum_absolute_energy_drift_hartree) <= thresholds["maximum_gauge_error"],
        "four_convergence_levels_are_present": len(convergence) == 4,
        "convergence_timesteps_halve": all(abs(right / left - 0.5) <= 1.0e-15 for left, right in zip(V250_CONVERGENCE_DT_AU[:-1], V250_CONVERGENCE_DT_AU[1:])),
        "first_state_change_ratio_is_second_order": thresholds["minimum_second_order_ratio"] <= state_ratios[0] <= thresholds["maximum_second_order_ratio"],
        "second_state_change_ratio_is_second_order": thresholds["minimum_second_order_ratio"] <= state_ratios[1] <= thresholds["maximum_second_order_ratio"],
        "first_energy_drift_ratio_is_second_order": thresholds["minimum_second_order_ratio"] <= energy_ratios[0] <= thresholds["maximum_second_order_ratio"],
        "second_energy_drift_ratio_is_second_order": thresholds["minimum_second_order_ratio"] <= energy_ratios[1] <= thresholds["maximum_second_order_ratio"],
        "all_raw_overlaps_are_physical_contractions": all(step.transport_metrics["physically_consistent"] for step in all_canonical_steps),
        "all_raw_overlaps_are_trajectory_quality": all(step.transport_metrics["trajectory_ready"] for step in all_canonical_steps),
        "all_polar_transports_are_unitary": max(transport_unitarity) <= thresholds["matrix_residual"],
        "all_step_singular_values_are_policy_retained": min(float(np.min(step.singular_values)) for step in all_canonical_steps) >= thresholds["minimum_contractive_singular_value"],
        "contractive_overlap_is_retained_as_evidence": all(np.allclose(step.overlap_start_end, 0.97 * np.eye(4)) for step in contraction.steps),
        "contractive_overlap_polar_factor_is_identity": all(np.allclose(step.transport_end_to_start, np.eye(4)) for step in contraction.steps),
        "contractive_overlap_does_not_change_norm": contraction.maximum_norm_drift <= thresholds["maximum_norm_drift"],
        "zero_soc_coordinate_equivalence": np.array_equal(zero_enabled.final_state.q, zero_disabled.final_state.q),
        "zero_soc_momentum_equivalence": np.array_equal(zero_enabled.final_state.p, zero_disabled.final_state.p),
        "zero_soc_spinor_equivalence": np.array_equal(zero_enabled.final_state.electronic_coefficients, zero_disabled.final_state.electronic_coefficients),
        "restricted_tdvp_claim_is_true": V250_TRAJECTORY_CLAIMS["restricted_single_packet_tdvp_validated"] is True,
        "full_multi_gaussian_tdvp_claim_is_false": V250_TRAJECTORY_CLAIMS["full_multi_gaussian_tdvp_validated"] is False,
        "real_pyscf_trajectory_claim_is_false": V250_TRAJECTORY_CLAIMS["real_pyscf_soc_trajectory_admitted"] is False,
        "general_accuracy_claim_is_false": V250_TRAJECTORY_CLAIMS["general_ab_initio_soc_dynamics_accuracy_validated"] is False,
    }
    checks = {name: bool(value) for name, value in checks.items()}
    if len(checks) != 45:
        raise AssertionError("v0.25.0 validation must define exactly 45 gates.")
    metrics = {
        "even_maximum_norm_drift": even.maximum_norm_drift,
        "even_maximum_energy_drift_hartree": even.maximum_absolute_energy_drift_hartree,
        "odd_maximum_norm_drift": odd.maximum_norm_drift,
        "odd_maximum_energy_drift_hartree": odd.maximum_absolute_energy_drift_hartree,
        "even_reversibility": {
            "q": float(np.max(np.abs(even_reverse.final_state.q - initial.q))),
            "p": float(np.max(np.abs(even_reverse.final_state.p - initial.p))),
            "spinor": float(np.max(np.abs(even_reverse.final_state.electronic_coefficients - initial.electronic_coefficients))),
        },
        "odd_reversibility": {
            "q": float(np.max(np.abs(odd_reverse.final_state.q - initial.q))),
            "p": float(np.max(np.abs(odd_reverse.final_state.p - initial.p))),
            "spinor": float(np.max(np.abs(odd_reverse.final_state.electronic_coefficients - initial.electronic_coefficients))),
        },
        "gauge_errors": {
            "q": float(np.max(np.abs(gauge_trajectory.final_state.q - odd.final_state.q))),
            "p": float(np.max(np.abs(gauge_trajectory.final_state.p - odd.final_state.p))),
            "spinor": gauge_coefficient_error,
        },
        "convergence_state_changes": state_changes,
        "convergence_state_change_ratios": state_ratios,
        "convergence_energy_drifts_hartree": energy_drifts,
        "convergence_energy_drift_ratios": energy_ratios,
        "maximum_transport_unitarity_residual": max(transport_unitarity),
        "minimum_retained_singular_value": min(float(np.min(step.singular_values)) for step in all_canonical_steps),
        "gate_count": len(checks),
    }
    audit = VariationalSOCValidationAuditV250(
        checks=checks,
        metrics=metrics,
        thresholds=thresholds,
        passed=bool(all(checks.values())),
    ).validate()
    return VariationalSOCValidationEvidenceV250(
        even_trajectory=even,
        even_reverse=even_reverse,
        odd_trajectory=odd,
        odd_reverse=odd_reverse,
        gauge_trajectory=gauge_trajectory,
        contraction_trajectory=contraction,
        zero_soc_enabled_trajectory=zero_enabled,
        zero_soc_disabled_trajectory=zero_disabled,
        convergence_trajectories=convergence,
        audit=audit,
    ).validate()


def save_variational_soc_validation_evidence_v250(path, evidence=None):
    from .campaign_io import save_campaign_json

    evidence = (
        run_variational_soc_validation_evidence_v250()
        if evidence is None
        else evidence.validate()
    )
    return save_campaign_json(path, evidence.as_dict())
