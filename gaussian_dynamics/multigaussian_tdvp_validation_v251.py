"""Deterministic validation evidence for the v0.25.1 TDVP metric layer."""

from dataclasses import dataclass

import numpy as np

from .analytic_soc_models_v220 import (
    AnalyticDoubletSOCProviderV220,
    AnalyticSingletTripletSOCProviderV220,
    DoubletSOCConfigV220,
)
from .complex_gauge_v21 import random_unitary_v21
from .multigaussian_tdvp_v251 import (
    MULTIGAUSSIAN_TDVP_ANSATZ_V251,
    MULTIGAUSSIAN_TDVP_SCHEMA_V251,
    VARIATIONAL_INTEGRATOR_V251,
    VARIATIONAL_METRIC_SOLVER_V251,
    VARIATIONAL_PRINCIPLE_V251,
    V251_TDVP_CLAIMS,
    FrozenGaussianSpinorStateV251,
    QuadraticSpinHamiltonianV251,
    _sha256_v251,
    build_variational_metric_system_v251,
    quadratic_spin_hamiltonian_from_provider_v251,
    reverse_frozen_width_multigaussian_tdvp_v251,
    run_frozen_width_multigaussian_tdvp_v251,
)


MULTIGAUSSIAN_TDVP_VALIDATION_SCHEMA_V251 = (
    "gnd-frozen-width-multigaussian-tdvp-validation-v0.25.1"
)
V251_CONVERGENCE_DT_AU = (0.1, 0.05, 0.025, 0.0125)
V251_CONVERGENCE_FINAL_TIME_AU = 0.4


def _initial_state_v251():
    coefficients = np.asarray(
        [
            [0.65 + 0.10j, 0.15 - 0.20j, 0.25 + 0.08j, -0.05j],
            [0.18 - 0.04j, -0.11 + 0.09j, 0.22 - 0.06j, 0.07 + 0.03j],
        ],
        dtype=complex,
    )
    return FrozenGaussianSpinorStateV251(
        q=np.asarray([-0.65, 0.75]),
        p=np.asarray([5.0, -3.0]),
        widths=np.asarray([2.6, 2.1]),
        coefficients=coefficients,
    ).normalized()


def _maximum_state_error_v251(left, right):
    return max(
        float(np.max(np.abs(left.q - right.q))),
        float(np.max(np.abs(left.p - right.p))),
        float(np.max(np.abs(left.widths - right.widths))),
        float(np.max(np.abs(left.coefficients - right.coefficients))),
        abs(float(left.time_au) - float(right.time_au)),
    )


def _trajectory_finite_v251(trajectory):
    states = [trajectory.initial_state] + [step.end for step in trajectory.steps]
    return bool(
        all(
            np.all(np.isfinite(state.q))
            and np.all(np.isfinite(state.p))
            and np.all(np.isfinite(state.coefficients))
            and np.isfinite(state.generalized_norm)
            for state in states
        )
    )


def _trajectory_receipt_v251(trajectory):
    return {
        "dt_au": (
            float(trajectory.steps[0].dt_au) if trajectory.steps else 0.0
        ),
        "step_count": len(trajectory.steps),
        "final_state": trajectory.final_state.as_dict(),
        "maximum_norm_drift": trajectory.maximum_norm_drift,
        "maximum_absolute_energy_drift_hartree": (
            trajectory.maximum_absolute_energy_drift_hartree
        ),
        "minimum_metric_rank": trajectory.minimum_metric_rank,
        "maximum_metric_nullity": trajectory.maximum_metric_nullity,
        "maximum_nonlinear_residual": trajectory.maximum_nonlinear_residual,
        "fingerprint": trajectory.fingerprint(),
    }


@dataclass(frozen=True)
class MultiGaussianTDVPValidationAuditV251:
    checks: dict
    metrics: dict
    thresholds: dict
    passed: bool

    def validate(self):
        if not isinstance(self.checks, dict) or len(self.checks) != 55:
            raise ValueError("v0.25.1 multi-Gaussian TDVP audit requires 55 gates.")
        if any(type(value) is not bool for value in self.checks.values()):
            raise TypeError("every v0.25.1 validation gate must be Boolean.")
        if type(self.passed) is not bool or self.passed != bool(
            all(self.checks.values())
        ):
            raise ValueError("v0.25.1 validation audit result is inconsistent.")
        _sha256_v251(self.metrics)
        _sha256_v251(self.thresholds)
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
class MultiGaussianTDVPValidationEvidenceV251:
    even_trajectory: object
    even_reverse: object
    odd_trajectory: object
    odd_reverse: object
    permutation_trajectory: object
    gauge_trajectory: object
    zero_soc_enabled_trajectory: object
    zero_soc_disabled_trajectory: object
    convergence_trajectories: tuple
    compatible_null_system: object
    harmonic_system: object
    audit: MultiGaussianTDVPValidationAuditV251

    def validate(self):
        for trajectory in (
            self.even_trajectory,
            self.even_reverse,
            self.odd_trajectory,
            self.odd_reverse,
            self.permutation_trajectory,
            self.gauge_trajectory,
            self.zero_soc_enabled_trajectory,
            self.zero_soc_disabled_trajectory,
            *self.convergence_trajectories,
        ):
            trajectory.validate()
        self.compatible_null_system.validate()
        self.harmonic_system.validate()
        self.audit.validate()
        if not self.audit.passed:
            failed = ", ".join(
                name for name, passed in self.audit.checks.items() if not passed
            )
            raise ValueError("v0.25.1 multi-Gaussian evidence failed: " + failed)
        return self

    @property
    def claims(self):
        return dict(V251_TDVP_CLAIMS)

    def as_dict(self):
        self.validate()
        return {
            "schema": MULTIGAUSSIAN_TDVP_VALIDATION_SCHEMA_V251,
            "trajectory_schema": MULTIGAUSSIAN_TDVP_SCHEMA_V251,
            "decisions": {
                "validated_ansatz": MULTIGAUSSIAN_TDVP_ANSATZ_V251,
                "variational_principle": VARIATIONAL_PRINCIPLE_V251,
                "integrator": VARIATIONAL_INTEGRATOR_V251,
                "metric_solver": VARIATIONAL_METRIC_SOLVER_V251,
                "widths": "positive, packet-specific, and frozen",
                "electronic_frame": "coordinate-independent complete spinor frame",
                "potential": "Hermitian matrix polynomial through degree two",
            },
            "canonical_trajectories": {
                "even": self.even_trajectory.as_dict(),
                "even_reverse": self.even_reverse.as_dict(),
                "odd": self.odd_trajectory.as_dict(),
                "odd_reverse": self.odd_reverse.as_dict(),
                "odd_gaussian_permutation": self.permutation_trajectory.as_dict(),
                "odd_constant_electronic_gauge": self.gauge_trajectory.as_dict(),
                "zero_soc_enabled": self.zero_soc_enabled_trajectory.as_dict(),
                "zero_soc_disabled": self.zero_soc_disabled_trajectory.as_dict(),
            },
            "convergence_receipts": [
                _trajectory_receipt_v251(trajectory)
                for trajectory in self.convergence_trajectories
            ],
            "compatible_null_metric": self.compatible_null_system.as_dict(),
            "single_packet_harmonic_metric": self.harmonic_system.as_dict(),
            "audit": self.audit.as_dict(),
            "claims": self.claims,
        }

    def fingerprint(self):
        return _sha256_v251(self.as_dict())


def run_multigaussian_tdvp_validation_evidence_v251():
    """Run deterministic even/odd, covariance, null-space, and order evidence."""

    initial = _initial_state_v251()
    even_model = quadratic_spin_hamiltonian_from_provider_v251(
        AnalyticSingletTripletSOCProviderV220()
    )
    odd_model = quadratic_spin_hamiltonian_from_provider_v251(
        AnalyticDoubletSOCProviderV220()
    )
    even = run_frozen_width_multigaussian_tdvp_v251(
        initial, even_model, dt_au=0.1, steps=4
    )
    odd = run_frozen_width_multigaussian_tdvp_v251(
        initial, odd_model, dt_au=0.1, steps=4
    )
    even_reverse = reverse_frozen_width_multigaussian_tdvp_v251(even)
    odd_reverse = reverse_frozen_width_multigaussian_tdvp_v251(odd)

    order = np.asarray([1, 0], dtype=int)
    permutation = run_frozen_width_multigaussian_tdvp_v251(
        initial.permuted(order), odd_model, dt_au=0.1, steps=4
    )
    expected_permutation = odd.final_state.permuted(order)

    unitary = random_unitary_v21(4, 25101)
    gauge = run_frozen_width_multigaussian_tdvp_v251(
        initial.gauge_transformed(unitary),
        odd_model.gauge_transformed(unitary),
        dt_au=0.1,
        steps=4,
    )
    expected_gauge = odd.final_state.gauge_transformed(unitary)

    zero_enabled_model = quadratic_spin_hamiltonian_from_provider_v251(
        AnalyticDoubletSOCProviderV220(
            DoubletSOCConfigV220(soc_scale=0.0, soc_enabled=True)
        )
    )
    zero_disabled_model = quadratic_spin_hamiltonian_from_provider_v251(
        AnalyticDoubletSOCProviderV220(
            DoubletSOCConfigV220(soc_scale=0.0, soc_enabled=False)
        )
    )
    zero_enabled = run_frozen_width_multigaussian_tdvp_v251(
        initial, zero_enabled_model, dt_au=0.1, steps=2
    )
    zero_disabled = run_frozen_width_multigaussian_tdvp_v251(
        initial, zero_disabled_model, dt_au=0.1, steps=2
    )

    convergence = tuple(
        run_frozen_width_multigaussian_tdvp_v251(
            initial,
            odd_model,
            dt_au=dt,
            steps=int(round(V251_CONVERGENCE_FINAL_TIME_AU / dt)),
        )
        for dt in V251_CONVERGENCE_DT_AU
    )
    state_changes = [
        _maximum_state_error_v251(left.final_state, right.final_state)
        for left, right in zip(convergence[:-1], convergence[1:])
    ]
    state_ratios = [
        right / left for left, right in zip(state_changes[:-1], state_changes[1:])
    ]

    harmonic_mass = 900.0
    harmonic_omega = 0.003
    harmonic_force = harmonic_mass * harmonic_omega**2
    harmonic_model = QuadraticSpinHamiltonianV251(
        harmonic_mass,
        np.zeros((1, 1)),
        np.zeros((1, 1)),
        np.asarray([[0.5 * harmonic_force]]),
        label="single-packet harmonic reduction",
    ).validate()
    harmonic_state = FrozenGaussianSpinorStateV251(
        q=[-0.4],
        p=[2.3],
        widths=[harmonic_mass * harmonic_omega],
        coefficients=[[1.0 + 0.0j]],
    ).normalized()
    harmonic_system = build_variational_metric_system_v251(
        harmonic_state, harmonic_model
    )

    compatible_state = FrozenGaussianSpinorStateV251(
        q=[-0.2, -0.2],
        p=[0.7, 0.7],
        widths=[2.4, 2.4],
        coefficients=[[0.6 + 0.1j], [0.3 - 0.2j]],
    ).normalized()
    compatible_null = build_variational_metric_system_v251(
        compatible_state, harmonic_model
    )

    thresholds = {
        "maximum_norm_drift": 2.0e-12,
        "maximum_energy_drift_hartree": 2.0e-11,
        "maximum_reversibility_error": 3.0e-12,
        "maximum_covariance_error": 3.0e-12,
        "maximum_nonlinear_residual": 2.0e-12,
        "maximum_null_rhs_relative": 2.0e-12,
        "maximum_linear_residual_relative": 2.0e-12,
        "minimum_second_order_ratio": 0.245,
        "maximum_second_order_ratio": 0.255,
        "harmonic_reduction_error": 3.0e-14,
    }
    even_return = even_reverse.final_state
    odd_return = odd_reverse.final_state
    permutation_errors = {
        "q": float(np.max(np.abs(permutation.final_state.q - expected_permutation.q))),
        "p": float(np.max(np.abs(permutation.final_state.p - expected_permutation.p))),
        "coefficients": float(
            np.max(
                np.abs(
                    permutation.final_state.coefficients
                    - expected_permutation.coefficients
                )
            )
        ),
    }
    gauge_errors = {
        "q": float(np.max(np.abs(gauge.final_state.q - expected_gauge.q))),
        "p": float(np.max(np.abs(gauge.final_state.p - expected_gauge.p))),
        "coefficients": float(
            np.max(np.abs(gauge.final_state.coefficients - expected_gauge.coefficients))
        ),
    }
    zero_models_equal = bool(
        np.array_equal(zero_enabled_model.H0, zero_disabled_model.H0)
        and np.array_equal(zero_enabled_model.H1, zero_disabled_model.H1)
        and np.array_equal(zero_enabled_model.H2, zero_disabled_model.H2)
    )
    harmonic_q_error = abs(
        harmonic_system.velocity[-2] - harmonic_state.p[0] / harmonic_mass
    )
    harmonic_p_error = abs(
        harmonic_system.velocity[-1] + harmonic_force * harmonic_state.q[0]
    )

    checks = {
        "validation_schema_is_v0251": MULTIGAUSSIAN_TDVP_VALIDATION_SCHEMA_V251.endswith("v0.25.1"),
        "frozen_width_multigaussian_ansatz_is_explicit": "coupled multi-Gaussian" in MULTIGAUSSIAN_TDVP_ANSATZ_V251 and "fixed-width" in MULTIGAUSSIAN_TDVP_ANSATZ_V251,
        "real_parameter_mclachlan_principle_is_frozen": "McLachlan" in VARIATIONAL_PRINCIPLE_V251,
        "fully_implicit_midpoint_is_frozen": "fully implicit midpoint" in VARIATIONAL_INTEGRATOR_V251,
        "full_svd_pseudoinverse_is_frozen": "full SVD" in VARIATIONAL_METRIC_SOLVER_V251,
        "packet_widths_remain_frozen": all(np.array_equal(step.start.widths, step.end.widths) for step in (*even.steps, *odd.steps)),
        "nuclear_scope_is_one_dimensional": all(state.q.ndim == 1 for state in (even.initial_state, odd.initial_state)),
        "electronic_frame_is_fixed": even.model.source["fixed_frame_verified"] is True and odd.model.source["fixed_frame_verified"] is True,
        "even_trajectory_is_finite": _trajectory_finite_v251(even),
        "even_generalized_norm_is_preserved": even.maximum_norm_drift <= thresholds["maximum_norm_drift"],
        "even_energy_drift_is_bounded": even.maximum_absolute_energy_drift_hartree <= thresholds["maximum_energy_drift_hartree"],
        "even_metric_rank_is_full": even.minimum_metric_rank == initial.parameter_count,
        "even_nonlinear_residual_is_certified": even.maximum_nonlinear_residual <= thresholds["maximum_nonlinear_residual"],
        "even_reverse_coordinate_returns": np.max(np.abs(even_return.q - initial.q)) <= thresholds["maximum_reversibility_error"],
        "even_reverse_momentum_returns": np.max(np.abs(even_return.p - initial.p)) <= thresholds["maximum_reversibility_error"],
        "even_reverse_spinor_returns": np.max(np.abs(even_return.coefficients - initial.coefficients)) <= thresholds["maximum_reversibility_error"],
        "even_reverse_time_returns": abs(even_return.time_au - initial.time_au) <= thresholds["maximum_reversibility_error"],
        "odd_trajectory_is_finite": _trajectory_finite_v251(odd),
        "odd_generalized_norm_is_preserved": odd.maximum_norm_drift <= thresholds["maximum_norm_drift"],
        "odd_energy_drift_is_bounded": odd.maximum_absolute_energy_drift_hartree <= thresholds["maximum_energy_drift_hartree"],
        "odd_metric_rank_is_full": odd.minimum_metric_rank == initial.parameter_count,
        "odd_nonlinear_residual_is_certified": odd.maximum_nonlinear_residual <= thresholds["maximum_nonlinear_residual"],
        "odd_reverse_coordinate_returns": np.max(np.abs(odd_return.q - initial.q)) <= thresholds["maximum_reversibility_error"],
        "odd_reverse_momentum_returns": np.max(np.abs(odd_return.p - initial.p)) <= thresholds["maximum_reversibility_error"],
        "odd_reverse_spinor_returns": np.max(np.abs(odd_return.coefficients - initial.coefficients)) <= thresholds["maximum_reversibility_error"],
        "odd_reverse_time_returns": abs(odd_return.time_au - initial.time_au) <= thresholds["maximum_reversibility_error"],
        "gaussian_permutation_coordinate_covariance": permutation_errors["q"] <= thresholds["maximum_covariance_error"],
        "gaussian_permutation_momentum_covariance": permutation_errors["p"] <= thresholds["maximum_covariance_error"],
        "gaussian_permutation_spinor_covariance": permutation_errors["coefficients"] <= thresholds["maximum_covariance_error"],
        "gaussian_permutation_energy_covariance": abs(permutation.maximum_absolute_energy_drift_hartree - odd.maximum_absolute_energy_drift_hartree) <= thresholds["maximum_covariance_error"],
        "constant_gauge_coordinate_covariance": gauge_errors["q"] <= thresholds["maximum_covariance_error"],
        "constant_gauge_momentum_covariance": gauge_errors["p"] <= thresholds["maximum_covariance_error"],
        "constant_gauge_spinor_covariance": gauge_errors["coefficients"] <= thresholds["maximum_covariance_error"],
        "constant_gauge_norm_covariance": abs(gauge.maximum_norm_drift - odd.maximum_norm_drift) <= thresholds["maximum_covariance_error"],
        "constant_gauge_energy_covariance": abs(gauge.maximum_absolute_energy_drift_hartree - odd.maximum_absolute_energy_drift_hartree) <= thresholds["maximum_covariance_error"],
        "compatible_duplicate_packets_are_rank_deficient": compatible_null.solve_receipt.rank < compatible_state.parameter_count,
        "compatible_duplicate_packets_expose_nullity": compatible_null.solve_receipt.nullity > 0,
        "compatible_null_rhs_projection_is_bounded": compatible_null.solve_receipt.null_rhs_relative <= thresholds["maximum_null_rhs_relative"],
        "compatible_pseudoinverse_residual_is_bounded": compatible_null.solve_receipt.linear_residual_relative <= thresholds["maximum_linear_residual_relative"],
        "single_packet_harmonic_coordinate_velocity_reduces": harmonic_q_error <= thresholds["harmonic_reduction_error"],
        "single_packet_harmonic_momentum_velocity_reduces": harmonic_p_error <= thresholds["harmonic_reduction_error"],
        "four_convergence_levels_are_present": len(convergence) == 4,
        "convergence_timesteps_halve": all(abs(right / left - 0.5) <= 1.0e-15 for left, right in zip(V251_CONVERGENCE_DT_AU[:-1], V251_CONVERGENCE_DT_AU[1:])),
        "first_state_change_ratio_is_second_order": thresholds["minimum_second_order_ratio"] <= state_ratios[0] <= thresholds["maximum_second_order_ratio"],
        "second_state_change_ratio_is_second_order": thresholds["minimum_second_order_ratio"] <= state_ratios[1] <= thresholds["maximum_second_order_ratio"],
        "zero_soc_models_are_exactly_equal": zero_models_equal,
        "zero_soc_coordinate_equivalence": np.array_equal(zero_enabled.final_state.q, zero_disabled.final_state.q),
        "zero_soc_momentum_equivalence": np.array_equal(zero_enabled.final_state.p, zero_disabled.final_state.p),
        "zero_soc_spinor_equivalence": np.array_equal(zero_enabled.final_state.coefficients, zero_disabled.final_state.coefficients),
        "frozen_width_multigaussian_claim_is_true": V251_TDVP_CLAIMS["frozen_width_multigaussian_tdvp_metric_validated"] is True,
        "adaptive_width_claim_is_false": V251_TDVP_CLAIMS["adaptive_gaussian_width_tdvp_validated"] is False,
        "spawning_claim_is_false": V251_TDVP_CLAIMS["dynamic_spawning_validated"] is False,
        "coordinate_dependent_gauge_claim_is_false": V251_TDVP_CLAIMS["coordinate_dependent_electronic_gauge_covariance_validated"] is False,
        "real_pyscf_trajectory_claim_is_false": V251_TDVP_CLAIMS["real_pyscf_soc_trajectory_admitted"] is False,
        "general_accuracy_claim_is_false": V251_TDVP_CLAIMS["general_ab_initio_soc_dynamics_accuracy_validated"] is False,
    }
    checks = {name: bool(value) for name, value in checks.items()}
    if len(checks) != 55:
        raise AssertionError("v0.25.1 validation must define exactly 55 gates.")

    metrics = {
        "even_maximum_norm_drift": even.maximum_norm_drift,
        "even_maximum_energy_drift_hartree": even.maximum_absolute_energy_drift_hartree,
        "odd_maximum_norm_drift": odd.maximum_norm_drift,
        "odd_maximum_energy_drift_hartree": odd.maximum_absolute_energy_drift_hartree,
        "even_reversibility": {
            "q": float(np.max(np.abs(even_return.q - initial.q))),
            "p": float(np.max(np.abs(even_return.p - initial.p))),
            "coefficients": float(np.max(np.abs(even_return.coefficients - initial.coefficients))),
        },
        "odd_reversibility": {
            "q": float(np.max(np.abs(odd_return.q - initial.q))),
            "p": float(np.max(np.abs(odd_return.p - initial.p))),
            "coefficients": float(np.max(np.abs(odd_return.coefficients - initial.coefficients))),
        },
        "permutation_errors": permutation_errors,
        "constant_gauge_errors": gauge_errors,
        "compatible_null_space": compatible_null.solve_receipt.as_dict(),
        "harmonic_reduction_errors": {"qdot": harmonic_q_error, "pdot": harmonic_p_error},
        "convergence_state_changes": state_changes,
        "convergence_state_change_ratios": state_ratios,
        "maximum_even_nonlinear_residual": even.maximum_nonlinear_residual,
        "maximum_odd_nonlinear_residual": odd.maximum_nonlinear_residual,
        "gate_count": len(checks),
    }
    audit = MultiGaussianTDVPValidationAuditV251(
        checks=checks,
        metrics=metrics,
        thresholds=thresholds,
        passed=bool(all(checks.values())),
    ).validate()
    return MultiGaussianTDVPValidationEvidenceV251(
        even_trajectory=even,
        even_reverse=even_reverse,
        odd_trajectory=odd,
        odd_reverse=odd_reverse,
        permutation_trajectory=permutation,
        gauge_trajectory=gauge,
        zero_soc_enabled_trajectory=zero_enabled,
        zero_soc_disabled_trajectory=zero_disabled,
        convergence_trajectories=convergence,
        compatible_null_system=compatible_null,
        harmonic_system=harmonic_system,
        audit=audit,
    ).validate()


def save_multigaussian_tdvp_validation_evidence_v251(path, evidence=None):
    from .campaign_io import save_campaign_json

    evidence = (
        run_multigaussian_tdvp_validation_evidence_v251()
        if evidence is None
        else evidence.validate()
    )
    return save_campaign_json(path, evidence.as_dict())
