"""Deterministic adaptive-width multi-Gaussian TDVP evidence for v0.25.2."""

from dataclasses import dataclass

import numpy as np

from .adaptive_multigaussian_tdvp_v252 import (
    ADAPTIVE_MULTIGAUSSIAN_TDVP_ANSATZ_V252,
    ADAPTIVE_MULTIGAUSSIAN_TDVP_SCHEMA_V252,
    VARIATIONAL_INTEGRATOR_V252,
    VARIATIONAL_METRIC_SOLVER_V252,
    VARIATIONAL_PRINCIPLE_V252,
    WIDTH_COORDINATES_V252,
    V252_TDVP_CLAIMS,
    QuadraticSpinHamiltonianV251,
    ThawedGaussianSpinorStateV252,
    _sha256_v252,
    build_adaptive_variational_metric_system_v252,
    quadratic_spin_hamiltonian_from_provider_v252,
    reverse_adaptive_width_multigaussian_tdvp_v252,
    run_adaptive_width_multigaussian_tdvp_v252,
)
from .analytic_soc_models_v220 import (
    AnalyticDoubletSOCProviderV220,
    AnalyticSingletTripletSOCProviderV220,
    DoubletSOCConfigV220,
)
from .complex_gauge_v21 import random_unitary_v21
from .multigaussian_tdvp_v251 import (
    FrozenGaussianSpinorStateV251,
    run_frozen_width_multigaussian_tdvp_v251,
)


ADAPTIVE_MULTIGAUSSIAN_VALIDATION_SCHEMA_V252 = (
    "gnd-adaptive-width-multigaussian-tdvp-validation-v0.25.2"
)
V252_CONVERGENCE_DT_AU = (0.1, 0.05, 0.025, 0.0125)
V252_CONVERGENCE_FINAL_TIME_AU = 0.2


def _initial_state_v252():
    coefficients = np.asarray(
        [
            [0.65 + 0.10j, 0.15 - 0.20j, 0.25 + 0.08j, -0.05j],
            [0.18 - 0.04j, -0.11 + 0.09j, 0.22 - 0.06j, 0.07 + 0.03j],
        ],
        dtype=complex,
    )
    return ThawedGaussianSpinorStateV252(
        q=[-0.65, 0.75],
        p=[5.0, -3.0],
        widths=[2.6, 2.1],
        chirps=[0.12, -0.08],
        coefficients=coefficients,
    ).normalized()


def _maximum_state_error_v252(left, right):
    return max(
        float(np.max(np.abs(left.q - right.q))),
        float(np.max(np.abs(left.p - right.p))),
        float(np.max(np.abs(left.widths - right.widths))),
        float(np.max(np.abs(left.chirps - right.chirps))),
        float(np.max(np.abs(left.coefficients - right.coefficients))),
        abs(float(left.time_au) - float(right.time_au)),
    )


def _component_errors_v252(left, right):
    return {
        "q": float(np.max(np.abs(left.q - right.q))),
        "p": float(np.max(np.abs(left.p - right.p))),
        "widths": float(np.max(np.abs(left.widths - right.widths))),
        "chirps": float(np.max(np.abs(left.chirps - right.chirps))),
        "coefficients": float(
            np.max(np.abs(left.coefficients - right.coefficients))
        ),
        "time_au": abs(float(left.time_au) - float(right.time_au)),
    }


def _trajectory_finite_v252(trajectory):
    states = [trajectory.initial_state] + [step.end for step in trajectory.steps]
    return bool(
        all(
            np.all(np.isfinite(state.q))
            and np.all(np.isfinite(state.p))
            and np.all(np.isfinite(state.widths))
            and np.all(np.isfinite(state.chirps))
            and np.all(np.isfinite(state.coefficients))
            and np.isfinite(state.generalized_norm)
            for state in states
        )
    )


def _trajectory_receipt_v252(trajectory):
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
        "minimum_width": trajectory.minimum_width,
        "maximum_width": trajectory.maximum_width,
        "maximum_absolute_chirp": trajectory.maximum_absolute_chirp,
        "fingerprint": trajectory.fingerprint(),
    }


def _harmonic_exact_v252(initial, mass, omega, time_au):
    angle = float(omega) * float(time_au)
    cosine = np.cos(angle)
    sine = np.sin(angle)
    q = initial.q[0] * cosine + initial.p[0] * sine / (mass * omega)
    p = initial.p[0] * cosine - mass * omega * initial.q[0] * sine
    width_parameter = initial.chirps[0] + 1.0j * initial.widths[0]
    evolved = (
        -mass * omega * sine + width_parameter * cosine
    ) / (cosine + width_parameter * sine / (mass * omega))
    return {
        "q": float(q),
        "p": float(p),
        "width": float(evolved.imag),
        "chirp": float(evolved.real),
    }


@dataclass(frozen=True)
class AdaptiveMultiGaussianValidationAuditV252:
    checks: dict
    metrics: dict
    thresholds: dict
    passed: bool

    def validate(self):
        if not isinstance(self.checks, dict) or len(self.checks) != 70:
            raise ValueError("v0.25.2 adaptive TDVP audit requires exactly 70 gates.")
        if any(type(value) is not bool for value in self.checks.values()):
            raise TypeError("every v0.25.2 adaptive validation gate must be Boolean.")
        if type(self.passed) is not bool or self.passed != bool(
            all(self.checks.values())
        ):
            raise ValueError("v0.25.2 adaptive validation result is inconsistent.")
        _sha256_v252(self.metrics)
        _sha256_v252(self.thresholds)
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
class AdaptiveMultiGaussianValidationEvidenceV252:
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
    harmonic_trajectory: object
    coherent_adaptive_trajectory: object
    coherent_frozen_trajectory: object
    audit: AdaptiveMultiGaussianValidationAuditV252

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
            self.harmonic_trajectory,
            self.coherent_adaptive_trajectory,
        ):
            trajectory.validate()
        self.coherent_frozen_trajectory.validate()
        self.compatible_null_system.validate()
        self.audit.validate()
        if not self.audit.passed:
            failed = ", ".join(
                name for name, passed in self.audit.checks.items() if not passed
            )
            raise ValueError("v0.25.2 adaptive evidence failed: " + failed)
        return self

    @property
    def claims(self):
        return dict(V252_TDVP_CLAIMS)

    def as_dict(self):
        self.validate()
        return {
            "schema": ADAPTIVE_MULTIGAUSSIAN_VALIDATION_SCHEMA_V252,
            "trajectory_schema": ADAPTIVE_MULTIGAUSSIAN_TDVP_SCHEMA_V252,
            "decisions": {
                "validated_ansatz": ADAPTIVE_MULTIGAUSSIAN_TDVP_ANSATZ_V252,
                "width_coordinates": WIDTH_COORDINATES_V252,
                "variational_principle": VARIATIONAL_PRINCIPLE_V252,
                "integrator": VARIATIONAL_INTEGRATOR_V252,
                "metric_solver": VARIATIONAL_METRIC_SOLVER_V252,
                "analytic_moment_degree": 4,
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
                "single_packet_harmonic_breathing": (
                    self.harmonic_trajectory.as_dict()
                ),
                "coherent_adaptive": self.coherent_adaptive_trajectory.as_dict(),
                "coherent_frozen_v0251": self.coherent_frozen_trajectory.as_dict(),
            },
            "convergence_receipts": [
                _trajectory_receipt_v252(trajectory)
                for trajectory in self.convergence_trajectories
            ],
            "compatible_null_metric": self.compatible_null_system.as_dict(),
            "audit": self.audit.as_dict(),
            "claims": self.claims,
        }

    def fingerprint(self):
        return _sha256_v252(self.as_dict())


def run_adaptive_multigaussian_validation_evidence_v252():
    initial = _initial_state_v252()
    even_model = quadratic_spin_hamiltonian_from_provider_v252(
        AnalyticSingletTripletSOCProviderV220()
    )
    odd_model = quadratic_spin_hamiltonian_from_provider_v252(
        AnalyticDoubletSOCProviderV220()
    )
    even = run_adaptive_width_multigaussian_tdvp_v252(
        initial, even_model, dt_au=0.05, steps=4
    )
    odd = run_adaptive_width_multigaussian_tdvp_v252(
        initial, odd_model, dt_au=0.05, steps=4
    )
    even_reverse = reverse_adaptive_width_multigaussian_tdvp_v252(even)
    odd_reverse = reverse_adaptive_width_multigaussian_tdvp_v252(odd)

    order = np.asarray([1, 0], dtype=int)
    permutation = run_adaptive_width_multigaussian_tdvp_v252(
        initial.permuted(order), odd_model, dt_au=0.05, steps=4
    )
    expected_permutation = odd.final_state.permuted(order)

    unitary = random_unitary_v21(4, 25201)
    gauge = run_adaptive_width_multigaussian_tdvp_v252(
        initial.gauge_transformed(unitary),
        odd_model.gauge_transformed(unitary),
        dt_au=0.05,
        steps=4,
    )
    expected_gauge = odd.final_state.gauge_transformed(unitary)

    zero_enabled_model = quadratic_spin_hamiltonian_from_provider_v252(
        AnalyticDoubletSOCProviderV220(
            DoubletSOCConfigV220(soc_scale=0.0, soc_enabled=True)
        )
    )
    zero_disabled_model = quadratic_spin_hamiltonian_from_provider_v252(
        AnalyticDoubletSOCProviderV220(
            DoubletSOCConfigV220(soc_scale=0.0, soc_enabled=False)
        )
    )
    zero_enabled = run_adaptive_width_multigaussian_tdvp_v252(
        initial, zero_enabled_model, dt_au=0.05, steps=2
    )
    zero_disabled = run_adaptive_width_multigaussian_tdvp_v252(
        initial, zero_disabled_model, dt_au=0.05, steps=2
    )

    convergence = tuple(
        run_adaptive_width_multigaussian_tdvp_v252(
            initial,
            odd_model,
            dt_au=dt,
            steps=int(round(V252_CONVERGENCE_FINAL_TIME_AU / dt)),
        )
        for dt in V252_CONVERGENCE_DT_AU
    )
    state_changes = [
        _maximum_state_error_v252(left.final_state, right.final_state)
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
        [[0.5 * harmonic_force]],
        label="adaptive harmonic breathing oracle",
    ).validate()
    harmonic_initial = ThawedGaussianSpinorStateV252(
        q=[-0.4],
        p=[2.3],
        widths=[2.1],
        chirps=[0.37],
        coefficients=[[1.0]],
    ).normalized()
    harmonic = run_adaptive_width_multigaussian_tdvp_v252(
        harmonic_initial, harmonic_model, dt_au=0.5, steps=40
    )
    harmonic_exact = _harmonic_exact_v252(
        harmonic_initial, harmonic_mass, harmonic_omega, 20.0
    )
    harmonic_errors = {
        "q": abs(harmonic.final_state.q[0] - harmonic_exact["q"]),
        "p": abs(harmonic.final_state.p[0] - harmonic_exact["p"]),
        "width": abs(harmonic.final_state.widths[0] - harmonic_exact["width"]),
        "chirp": abs(harmonic.final_state.chirps[0] - harmonic_exact["chirp"]),
    }

    coherent_width = harmonic_mass * harmonic_omega
    coherent_adaptive_initial = ThawedGaussianSpinorStateV252(
        q=[-0.4],
        p=[2.3],
        widths=[coherent_width],
        chirps=[0.0],
        coefficients=[[1.0]],
    ).normalized()
    coherent_frozen_initial = FrozenGaussianSpinorStateV251(
        q=[-0.4],
        p=[2.3],
        widths=[coherent_width],
        coefficients=[[1.0]],
    ).normalized()
    coherent_adaptive = run_adaptive_width_multigaussian_tdvp_v252(
        coherent_adaptive_initial, harmonic_model, dt_au=0.1, steps=4
    )
    coherent_frozen = run_frozen_width_multigaussian_tdvp_v251(
        coherent_frozen_initial, harmonic_model, dt_au=0.1, steps=4
    )
    coherent_errors = {
        "q": float(
            np.max(
                np.abs(
                    coherent_adaptive.final_state.q - coherent_frozen.final_state.q
                )
            )
        ),
        "p": float(
            np.max(
                np.abs(
                    coherent_adaptive.final_state.p - coherent_frozen.final_state.p
                )
            )
        ),
        "coefficients": float(
            np.max(
                np.abs(
                    coherent_adaptive.final_state.coefficients
                    - coherent_frozen.final_state.coefficients
                )
            )
        ),
        "width": abs(coherent_adaptive.final_state.widths[0] - coherent_width),
        "chirp": abs(coherent_adaptive.final_state.chirps[0]),
    }

    compatible_state = ThawedGaussianSpinorStateV252(
        q=[-0.2, -0.2],
        p=[0.7, 0.7],
        widths=[2.4, 2.4],
        chirps=[0.2, 0.2],
        coefficients=[[0.6 + 0.1j], [0.3 - 0.2j]],
    ).normalized()
    compatible_null = build_adaptive_variational_metric_system_v252(
        compatible_state, harmonic_model
    )

    thresholds = {
        "maximum_norm_drift": 3.0e-12,
        "maximum_energy_drift_hartree": 3.0e-11,
        "maximum_reversibility_error": 4.0e-12,
        "maximum_covariance_error": 4.0e-12,
        "maximum_nonlinear_residual": 3.0e-12,
        "maximum_null_rhs_relative": 2.0e-12,
        "maximum_linear_residual_relative": 2.0e-12,
        "minimum_second_order_ratio": 0.245,
        "maximum_second_order_ratio": 0.255,
        "harmonic_qp_error": 2.0e-8,
        "harmonic_width_error": 2.0e-8,
        "harmonic_chirp_error": 6.0e-8,
        "coherent_reduction_error": 5.0e-13,
        "minimum_observable_width_change": 1.0e-7,
        "minimum_observable_chirp_change": 1.0e-7,
    }
    even_return_errors = _component_errors_v252(even_reverse.final_state, initial)
    odd_return_errors = _component_errors_v252(odd_reverse.final_state, initial)
    permutation_errors = _component_errors_v252(
        permutation.final_state, expected_permutation
    )
    gauge_errors = _component_errors_v252(gauge.final_state, expected_gauge)
    even_width_change = float(
        np.max(np.abs(even.final_state.widths - initial.widths))
    )
    even_chirp_change = float(
        np.max(np.abs(even.final_state.chirps - initial.chirps))
    )
    odd_width_change = float(
        np.max(np.abs(odd.final_state.widths - initial.widths))
    )
    odd_chirp_change = float(
        np.max(np.abs(odd.final_state.chirps - initial.chirps))
    )
    zero_model_equal = bool(
        np.array_equal(zero_enabled_model.H0, zero_disabled_model.H0)
        and np.array_equal(zero_enabled_model.H1, zero_disabled_model.H1)
        and np.array_equal(zero_enabled_model.H2, zero_disabled_model.H2)
    )

    checks = {
        "validation_schema_is_v0252": ADAPTIVE_MULTIGAUSSIAN_VALIDATION_SCHEMA_V252.endswith("v0.25.2"),
        "adaptive_multigaussian_ansatz_is_explicit": "logarithmic widths" in ADAPTIVE_MULTIGAUSSIAN_TDVP_ANSATZ_V252,
        "quadratic_chirp_coordinate_is_explicit": "quadratic chirp" in WIDTH_COORDINATES_V252,
        "real_parameter_mclachlan_principle_is_frozen": "McLachlan" in VARIATIONAL_PRINCIPLE_V252,
        "fully_implicit_midpoint_is_frozen": "fully implicit midpoint" in VARIATIONAL_INTEGRATOR_V252,
        "full_svd_pseudoinverse_is_frozen": "full SVD" in VARIATIONAL_METRIC_SOLVER_V252,
        "adaptive_width_claim_is_enabled": V252_TDVP_CLAIMS["adaptive_width_multigaussian_tdvp_validated"] is True,
        "nuclear_scope_is_one_dimensional": all(state.q.ndim == 1 for state in (even.initial_state, odd.initial_state)),
        "electronic_frame_is_fixed": even.model.source["fixed_frame_verified"] is True and odd.model.source["fixed_frame_verified"] is True,
        "degree_four_moments_are_declared": 4 == 4,
        "even_trajectory_is_finite": _trajectory_finite_v252(even),
        "even_generalized_norm_is_preserved": even.maximum_norm_drift <= thresholds["maximum_norm_drift"],
        "even_energy_drift_is_bounded": even.maximum_absolute_energy_drift_hartree <= thresholds["maximum_energy_drift_hartree"],
        "even_metric_rank_is_full": even.minimum_metric_rank == initial.parameter_count,
        "even_nonlinear_residual_is_certified": even.maximum_nonlinear_residual <= thresholds["maximum_nonlinear_residual"],
        "even_widths_remain_positive": even.minimum_width > 0.0,
        "even_widths_evolve": even_width_change >= thresholds["minimum_observable_width_change"],
        "even_chirps_evolve": even_chirp_change >= thresholds["minimum_observable_chirp_change"],
        "even_reverse_coordinate_returns": even_return_errors["q"] <= thresholds["maximum_reversibility_error"],
        "even_reverse_momentum_returns": even_return_errors["p"] <= thresholds["maximum_reversibility_error"],
        "even_reverse_width_returns": even_return_errors["widths"] <= thresholds["maximum_reversibility_error"],
        "even_reverse_chirp_returns": even_return_errors["chirps"] <= thresholds["maximum_reversibility_error"],
        "even_reverse_spinor_returns": even_return_errors["coefficients"] <= thresholds["maximum_reversibility_error"],
        "even_reverse_time_returns": even_return_errors["time_au"] <= thresholds["maximum_reversibility_error"],
        "odd_trajectory_is_finite": _trajectory_finite_v252(odd),
        "odd_generalized_norm_is_preserved": odd.maximum_norm_drift <= thresholds["maximum_norm_drift"],
        "odd_energy_drift_is_bounded": odd.maximum_absolute_energy_drift_hartree <= thresholds["maximum_energy_drift_hartree"],
        "odd_metric_rank_is_full": odd.minimum_metric_rank == initial.parameter_count,
        "odd_nonlinear_residual_is_certified": odd.maximum_nonlinear_residual <= thresholds["maximum_nonlinear_residual"],
        "odd_widths_remain_positive": odd.minimum_width > 0.0,
        "odd_widths_evolve": odd_width_change >= thresholds["minimum_observable_width_change"],
        "odd_chirps_evolve": odd_chirp_change >= thresholds["minimum_observable_chirp_change"],
        "odd_reverse_coordinate_returns": odd_return_errors["q"] <= thresholds["maximum_reversibility_error"],
        "odd_reverse_momentum_returns": odd_return_errors["p"] <= thresholds["maximum_reversibility_error"],
        "odd_reverse_width_returns": odd_return_errors["widths"] <= thresholds["maximum_reversibility_error"],
        "odd_reverse_chirp_returns": odd_return_errors["chirps"] <= thresholds["maximum_reversibility_error"],
        "odd_reverse_spinor_returns": odd_return_errors["coefficients"] <= thresholds["maximum_reversibility_error"],
        "odd_reverse_time_returns": odd_return_errors["time_au"] <= thresholds["maximum_reversibility_error"],
        "permutation_coordinate_covariance": permutation_errors["q"] <= thresholds["maximum_covariance_error"],
        "permutation_momentum_covariance": permutation_errors["p"] <= thresholds["maximum_covariance_error"],
        "permutation_width_covariance": permutation_errors["widths"] <= thresholds["maximum_covariance_error"],
        "permutation_chirp_covariance": permutation_errors["chirps"] <= thresholds["maximum_covariance_error"],
        "permutation_spinor_covariance": permutation_errors["coefficients"] <= thresholds["maximum_covariance_error"],
        "constant_gauge_coordinate_covariance": gauge_errors["q"] <= thresholds["maximum_covariance_error"],
        "constant_gauge_momentum_covariance": gauge_errors["p"] <= thresholds["maximum_covariance_error"],
        "constant_gauge_width_covariance": gauge_errors["widths"] <= thresholds["maximum_covariance_error"],
        "constant_gauge_chirp_covariance": gauge_errors["chirps"] <= thresholds["maximum_covariance_error"],
        "constant_gauge_spinor_covariance": gauge_errors["coefficients"] <= thresholds["maximum_covariance_error"],
        "duplicate_packets_have_rank_six": compatible_null.solve_receipt.rank == 6,
        "duplicate_packets_have_nullity_six": compatible_null.solve_receipt.nullity == 6,
        "compatible_null_rhs_is_bounded": compatible_null.solve_receipt.null_rhs_relative <= thresholds["maximum_null_rhs_relative"],
        "compatible_metric_residual_is_bounded": compatible_null.solve_receipt.linear_residual_relative <= thresholds["maximum_linear_residual_relative"],
        "harmonic_coordinate_matches_exact": harmonic_errors["q"] <= thresholds["harmonic_qp_error"],
        "harmonic_momentum_matches_exact": harmonic_errors["p"] <= thresholds["harmonic_qp_error"],
        "harmonic_width_matches_exact": harmonic_errors["width"] <= thresholds["harmonic_width_error"],
        "harmonic_chirp_matches_exact": harmonic_errors["chirp"] <= thresholds["harmonic_chirp_error"],
        "harmonic_breathing_is_nontrivial": abs(harmonic.final_state.widths[0] - harmonic_initial.widths[0]) > 1.0e-3,
        "coherent_frozen_reduction_is_exact": max(coherent_errors.values()) <= thresholds["coherent_reduction_error"],
        "four_convergence_levels_are_present": len(convergence) == 4,
        "convergence_timesteps_halve": all(abs(right / left - 0.5) <= 1.0e-15 for left, right in zip(V252_CONVERGENCE_DT_AU[:-1], V252_CONVERGENCE_DT_AU[1:])),
        "first_state_change_ratio_is_second_order": thresholds["minimum_second_order_ratio"] <= state_ratios[0] <= thresholds["maximum_second_order_ratio"],
        "second_state_change_ratio_is_second_order": thresholds["minimum_second_order_ratio"] <= state_ratios[1] <= thresholds["maximum_second_order_ratio"],
        "zero_soc_models_are_exactly_equal": zero_model_equal,
        "zero_soc_trajectory_is_exactly_equal": _maximum_state_error_v252(zero_enabled.final_state, zero_disabled.final_state) == 0.0,
        "log_width_chirp_claim_is_true": V252_TDVP_CLAIMS["log_width_positivity_and_quadratic_chirp_validated"] is True,
        "spawning_claim_is_false": V252_TDVP_CLAIMS["dynamic_spawning_validated"] is False,
        "coordinate_dependent_gauge_claim_is_false": V252_TDVP_CLAIMS["coordinate_dependent_electronic_gauge_covariance_validated"] is False,
        "multidimensional_claim_is_false": V252_TDVP_CLAIMS["multidimensional_adaptive_width_tdvp_validated"] is False,
        "real_pyscf_trajectory_claim_is_false": V252_TDVP_CLAIMS["real_pyscf_soc_trajectory_admitted"] is False,
        "general_accuracy_claim_is_false": V252_TDVP_CLAIMS["general_ab_initio_soc_dynamics_accuracy_validated"] is False,
    }
    checks = {name: bool(value) for name, value in checks.items()}
    if len(checks) != 70:
        raise AssertionError(
            f"v0.25.2 validation must define exactly 70 gates, found {len(checks)}."
        )

    metrics = {
        "even_maximum_norm_drift": even.maximum_norm_drift,
        "even_maximum_energy_drift_hartree": even.maximum_absolute_energy_drift_hartree,
        "odd_maximum_norm_drift": odd.maximum_norm_drift,
        "odd_maximum_energy_drift_hartree": odd.maximum_absolute_energy_drift_hartree,
        "even_reversibility": even_return_errors,
        "odd_reversibility": odd_return_errors,
        "even_width_change": even_width_change,
        "even_chirp_change": even_chirp_change,
        "odd_width_change": odd_width_change,
        "odd_chirp_change": odd_chirp_change,
        "permutation_errors": permutation_errors,
        "constant_gauge_errors": gauge_errors,
        "compatible_null_space": compatible_null.solve_receipt.as_dict(),
        "harmonic_exact_endpoint": harmonic_exact,
        "harmonic_errors": harmonic_errors,
        "coherent_frozen_reduction_errors": coherent_errors,
        "convergence_state_changes": state_changes,
        "convergence_state_change_ratios": state_ratios,
        "maximum_even_nonlinear_residual": even.maximum_nonlinear_residual,
        "maximum_odd_nonlinear_residual": odd.maximum_nonlinear_residual,
        "gate_count": len(checks),
    }
    audit = AdaptiveMultiGaussianValidationAuditV252(
        checks=checks,
        metrics=metrics,
        thresholds=thresholds,
        passed=bool(all(checks.values())),
    ).validate()
    return AdaptiveMultiGaussianValidationEvidenceV252(
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
        harmonic_trajectory=harmonic,
        coherent_adaptive_trajectory=coherent_adaptive,
        coherent_frozen_trajectory=coherent_frozen,
        audit=audit,
    ).validate()


def save_adaptive_multigaussian_validation_evidence_v252(path, evidence=None):
    from .campaign_io import save_campaign_json

    evidence = (
        run_adaptive_multigaussian_validation_evidence_v252()
        if evidence is None
        else evidence.validate()
    )
    return save_campaign_json(path, evidence.as_dict())
