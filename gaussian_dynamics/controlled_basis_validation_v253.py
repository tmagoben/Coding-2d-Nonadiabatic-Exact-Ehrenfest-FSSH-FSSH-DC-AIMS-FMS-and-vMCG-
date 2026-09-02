"""Deterministic scientific evidence for controlled basis adaptation v0.25.3."""

from dataclasses import dataclass

import numpy as np

from .adaptive_multigaussian_tdvp_v252 import (
    QuadraticSpinHamiltonianV252,
    ThawedGaussianSpinorStateV252,
    _sha256_v252,
    adaptive_implicit_midpoint_tdvp_step_v252,
    pack_adaptive_variational_parameters_v252,
    quadratic_spin_hamiltonian_from_provider_v252,
)
from .analytic_soc_models_v220 import (
    AnalyticDoubletSOCProviderV220,
    AnalyticSingletTripletSOCProviderV220,
    DoubletSOCConfigV220,
)
from .complex_gauge_v21 import random_unitary_v21
from .controlled_basis_adaptation_v253 import (
    CONTROLLED_BASIS_SCHEMA_V253,
    EVENT_ORDER_V253,
    POTENTIAL_CONTRACT_V253,
    PROJECTION_POLICY_V253,
    SPAWN_SCORE_V253,
    V253_CONTROLLED_BASIS_CLAIMS,
    ControlledBasisSettingsV253,
    SpawnCandidateV253,
    adapt_basis_once_v253,
    coefficient_activation_implicit_step_v253,
    evaluate_spawn_candidate_v253,
    generate_spawn_candidates_v253,
    run_controlled_basis_dynamics_v253,
)


CONTROLLED_BASIS_VALIDATION_SCHEMA_V253 = (
    "gnd-controlled-basis-validation-v0.25.3"
)


def _initial_state_v253():
    return ThawedGaussianSpinorStateV252(
        q=[-0.65, 0.75],
        p=[5.0, -3.0],
        widths=[2.6, 2.1],
        chirps=[0.12, -0.08],
        coefficients=[
            [0.65 + 0.10j, 0.15 - 0.20j, 0.25 + 0.08j, -0.05j],
            [0.18 - 0.04j, -0.11 + 0.09j, 0.22 - 0.06j, 0.07 + 0.03j],
        ],
    ).normalized()


def _parameter_error_v253(left, right):
    if left.ngaussian != right.ngaussian or left.nstate != right.nstate:
        return float("inf")
    return float(
        np.max(
            np.abs(
                pack_adaptive_variational_parameters_v252(left)
                - pack_adaptive_variational_parameters_v252(right)
            )
        )
    )


@dataclass(frozen=True)
class ControlledBasisValidationAuditV253:
    checks: dict
    metrics: dict
    thresholds: dict
    passed: bool

    def validate(self):
        if not isinstance(self.checks, dict) or len(self.checks) != 60:
            raise ValueError("v0.25.3 controlled-basis audit requires exactly 60 gates.")
        if any(type(value) is not bool for value in self.checks.values()):
            raise TypeError("every v0.25.3 controlled-basis gate must be Boolean.")
        if type(self.passed) is not bool or self.passed != bool(all(self.checks.values())):
            raise ValueError("v0.25.3 controlled-basis audit result is inconsistent.")
        _sha256_v252(self.metrics)
        _sha256_v252(self.thresholds)
        return self

    def as_dict(self):
        self.validate()
        return {
            "checks": dict(self.checks),
            "metrics": self.metrics,
            "thresholds": self.thresholds,
            "passed": self.passed,
        }


@dataclass(frozen=True)
class ControlledBasisValidationEvidenceV253:
    odd_spawn_event: object
    even_spawn_event: object
    permutation_spawn_event: object
    gauge_spawn_event: object
    activation_step: object
    prune_event: object
    merge_event: object
    controlled_trajectory: object
    no_event_trajectory: object
    no_event_reference: object
    maximum_count_event: object
    duplicate_candidate_evaluation: object
    zero_soc_enabled_evaluation: object
    zero_soc_disabled_evaluation: object
    audit: ControlledBasisValidationAuditV253

    @property
    def claims(self):
        return dict(V253_CONTROLLED_BASIS_CLAIMS)

    def validate(self):
        for event in (
            self.odd_spawn_event,
            self.even_spawn_event,
            self.permutation_spawn_event,
            self.gauge_spawn_event,
            self.prune_event,
            self.merge_event,
            self.maximum_count_event,
        ):
            event.validate()
        self.activation_step.validate()
        self.controlled_trajectory.validate()
        self.no_event_trajectory.validate()
        self.no_event_reference.validate()
        self.duplicate_candidate_evaluation.validate()
        self.zero_soc_enabled_evaluation.validate()
        self.zero_soc_disabled_evaluation.validate()
        self.audit.validate()
        if not self.audit.passed:
            failed = ", ".join(name for name, value in self.audit.checks.items() if not value)
            raise ValueError("v0.25.3 controlled-basis evidence failed: " + failed)
        return self

    def as_dict(self):
        self.validate()
        return {
            "schema": CONTROLLED_BASIS_VALIDATION_SCHEMA_V253,
            "trajectory_schema": CONTROLLED_BASIS_SCHEMA_V253,
            "decisions": {
                "residual_score": SPAWN_SCORE_V253,
                "projection_policy": PROJECTION_POLICY_V253,
                "event_order": EVENT_ORDER_V253,
                "potential_contract": POTENTIAL_CONTRACT_V253,
                "newborn_activation": (
                    "all electronic coefficients active immediately; q,p,log-width,chirp "
                    "frozen until coefficient row population reaches the gate"
                ),
            },
            "events": {
                "odd_spawn": self.odd_spawn_event.as_dict(),
                "even_spawn": self.even_spawn_event.as_dict(),
                "permutation_spawn": self.permutation_spawn_event.as_dict(),
                "constant_gauge_spawn": self.gauge_spawn_event.as_dict(),
                "prune": self.prune_event.as_dict(),
                "merge": self.merge_event.as_dict(),
                "maximum_count_noop": self.maximum_count_event.as_dict(),
            },
            "activation_step": self.activation_step.as_dict(),
            "controlled_trajectory": self.controlled_trajectory.as_dict(),
            "no_event_reduction": {
                "controlled": self.no_event_trajectory.as_dict(),
                "v0252_reference": self.no_event_reference.as_dict(),
            },
            "candidate_controls": {
                "duplicate": self.duplicate_candidate_evaluation.as_dict(),
                "zero_soc_enabled": self.zero_soc_enabled_evaluation.as_dict(),
                "zero_soc_disabled": self.zero_soc_disabled_evaluation.as_dict(),
            },
            "audit": self.audit.as_dict(),
            "claims": self.claims,
        }

    def fingerprint(self):
        return _sha256_v252(self.as_dict())


def run_controlled_basis_validation_evidence_v253():
    state = _initial_state_v253()
    odd_model = quadratic_spin_hamiltonian_from_provider_v252(
        AnalyticDoubletSOCProviderV220()
    )
    even_model = quadratic_spin_hamiltonian_from_provider_v252(
        AnalyticSingletTripletSOCProviderV220()
    )
    metadata = {
        "packet_ids": ("g000000", "g000001"),
        "packet_ages": (2, 2),
        "next_packet_serial": 2,
    }
    odd_spawn = adapt_basis_once_v253(state, odd_model, **metadata)
    even_spawn = adapt_basis_once_v253(state, even_model, **metadata)

    order = np.asarray([1, 0])
    permutation_spawn = adapt_basis_once_v253(
        state.permuted(order),
        odd_model,
        packet_ids=("g000001", "g000000"),
        packet_ages=(2, 2),
        next_packet_serial=2,
    )
    restored_permutation = permutation_spawn.after.permuted(np.asarray([1, 0, 2]))

    unitary = random_unitary_v21(4, 25301)
    gauge_spawn = adapt_basis_once_v253(
        state.gauge_transformed(unitary),
        odd_model.gauge_transformed(unitary),
        **metadata,
    )
    expected_gauge = odd_spawn.after.gauge_transformed(unitary)

    activation = coefficient_activation_implicit_step_v253(
        odd_spawn.after,
        odd_model,
        0.02,
        np.asarray([True, True, False]),
    )
    newborn_population = float(np.sum(np.abs(activation.end.coefficients[-1]) ** 2))
    dormant_shape_drift = max(
        abs(activation.end.q[-1] - activation.start.q[-1]),
        abs(activation.end.p[-1] - activation.start.p[-1]),
        abs(activation.end.widths[-1] - activation.start.widths[-1]),
        abs(activation.end.chirps[-1] - activation.start.chirps[-1]),
    )

    scalar_model = QuadraticSpinHamiltonianV252(
        900.0,
        np.zeros((1, 1)),
        np.zeros((1, 1)),
        np.asarray([[0.004]]),
        label="v0.25.3 projection lifecycle oracle",
    ).validate()
    prune_state = ThawedGaussianSpinorStateV252(
        q=[-1.0, 1.0], p=[0.0, 0.0], widths=[2.0, 2.0],
        chirps=[0.0, 0.0], coefficients=[[1.0], [1.0e-6]],
    ).normalized()
    prune = adapt_basis_once_v253(
        prune_state,
        scalar_model,
        packet_ids=("g000000", "g000001"),
        packet_ages=(64, 64),
        next_packet_serial=2,
    )
    merge_state = ThawedGaussianSpinorStateV252(
        q=[0.0, 0.01], p=[0.0, 0.0], widths=[2.0, 2.0],
        chirps=[0.0, 0.0], coefficients=[[0.7], [0.3]],
    ).normalized()
    merge_settings = ControlledBasisSettingsV253(
        maximum_merge_projection_loss=1.0e-4,
        maximum_event_energy_jump_hartree=1.0e-4,
    )
    merge = adapt_basis_once_v253(
        merge_state,
        scalar_model,
        packet_ids=("g000000", "g000001"),
        packet_ages=(2, 2),
        next_packet_serial=2,
        settings=merge_settings,
    )

    controlled = run_controlled_basis_dynamics_v253(
        state, odd_model, dt_au=0.02, steps=3
    )
    no_event_settings = ControlledBasisSettingsV253(
        spawn_residual_capture_threshold=1.0e3
    )
    no_event = run_controlled_basis_dynamics_v253(
        state, odd_model, dt_au=0.02, steps=1, settings=no_event_settings
    )
    no_event_reference = adaptive_implicit_midpoint_tdvp_step_v252(
        state, odd_model, 0.02
    )

    maximum_count_settings = ControlledBasisSettingsV253(maximum_packet_count=2)
    maximum_count = adapt_basis_once_v253(
        state, odd_model, settings=maximum_count_settings, **metadata
    )
    duplicate = SpawnCandidateV253(
        state.q[0], state.p[0], state.widths[0], state.chirps[0], 0, "external"
    )
    duplicate_evaluation = evaluate_spawn_candidate_v253(
        state, odd_model, duplicate
    )

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
    zero_candidate = generate_spawn_candidates_v253(state)[0]
    zero_enabled = evaluate_spawn_candidate_v253(
        state, zero_enabled_model, zero_candidate
    )
    zero_disabled = evaluate_spawn_candidate_v253(
        state, zero_disabled_model, zero_candidate
    )

    permutation_error = _parameter_error_v253(odd_spawn.after, restored_permutation)
    gauge_error = _parameter_error_v253(gauge_spawn.after, expected_gauge)
    no_event_error = _parameter_error_v253(
        no_event.final_state, no_event_reference.end
    )
    zero_score_error = abs(
        zero_enabled.residual_capture - zero_disabled.residual_capture
    )
    zero_coupling_error = float(
        np.max(
            np.abs(
                zero_enabled.orthogonalized_residual_coupling
                - zero_disabled.orthogonalized_residual_coupling
            )
        )
    )
    merge_overlap = float(abs(merge_state.nuclear_overlap_matrix()[0, 1]))
    thresholds = {
        "projection_loss": 2.0e-10,
        "projection_energy_jump_hartree": 2.0e-10,
        "covariance_error": 2.0e-11,
        "no_event_reduction_error": 5.0e-14,
        "activation_norm_drift": 3.0e-8,
        "activation_residual": 2.0e-10,
        "minimum_newborn_population": 1.0e-12,
        "dormant_shape_drift": 2.0e-14,
        "prune_loss": 2.0e-12,
        "prune_energy_jump_hartree": 2.0e-8,
        "merge_loss": 1.0e-4,
        "merge_energy_jump_hartree": 1.0e-4,
        "zero_soc_error": 2.0e-13,
    }
    metrics = {
        "odd_best_residual_capture": max(item.residual_capture for item in odd_spawn.candidate_evaluations),
        "even_best_residual_capture": max(item.residual_capture for item in even_spawn.candidate_evaluations),
        "odd_spawn_projection_loss": odd_spawn.projection.relative_projection_loss,
        "odd_spawn_energy_jump_hartree": odd_spawn.projection.energy_jump_hartree,
        "even_spawn_projection_loss": even_spawn.projection.relative_projection_loss,
        "even_spawn_energy_jump_hartree": even_spawn.projection.energy_jump_hartree,
        "permutation_error": permutation_error,
        "gauge_error": gauge_error,
        "newborn_population_after_one_activation_step": newborn_population,
        "dormant_shape_drift": dormant_shape_drift,
        "activation_norm_drift": abs(activation.norm_change),
        "activation_nonlinear_residual": activation.nonlinear_residual_norm,
        "prune_projection_loss": prune.projection.relative_projection_loss,
        "prune_energy_jump_hartree": prune.projection.energy_jump_hartree,
        "merge_pair_overlap": merge_overlap,
        "merge_projection_loss": merge.projection.relative_projection_loss,
        "merge_energy_jump_hartree": merge.projection.energy_jump_hartree,
        "controlled_maximum_norm_drift": controlled.maximum_norm_drift,
        "controlled_maximum_projection_loss": controlled.maximum_projection_loss,
        "no_event_reduction_error": no_event_error,
        "zero_soc_score_error": zero_score_error,
        "zero_soc_coupling_error": zero_coupling_error,
    }
    checks = {
        "validation_schema_is_v0253": CONTROLLED_BASIS_VALIDATION_SCHEMA_V253.endswith("v0.25.3"),
        "trajectory_schema_is_v0253": CONTROLLED_BASIS_SCHEMA_V253.endswith("v0.25.3"),
        "analytic_residual_score_is_frozen": "dPsi/dt+iHPsi" in SPAWN_SCORE_V253,
        "full_svd_projection_is_frozen": "full-SVD" in PROJECTION_POLICY_V253,
        "one_event_order_is_frozen": "at most one event" in EVENT_ORDER_V253,
        "quadratic_fixed_frame_scope_is_frozen": "quadratic" in POTENTIAL_CONTRACT_V253,
        "odd_complete_spin_model_admitted": odd_model.nstate == 4 and odd_model.complete_spin_manifold,
        "even_complete_spin_model_admitted": even_model.nstate == 4 and even_model.complete_spin_manifold,
        "odd_spawn_event_accepted": odd_spawn.event_kind == "spawn",
        "even_spawn_event_accepted": even_spawn.event_kind == "spawn",
        "odd_spawn_added_one_packet": odd_spawn.after.ngaussian == state.ngaussian + 1,
        "even_spawn_added_one_packet": even_spawn.after.ngaussian == state.ngaussian + 1,
        "odd_spawn_has_admitted_candidate": any(item.admitted for item in odd_spawn.candidate_evaluations),
        "even_spawn_has_admitted_candidate": any(item.admitted for item in even_spawn.candidate_evaluations),
        "odd_selected_candidate_is_highest_score": odd_spawn.selected_candidate.canonical_key() == sorted([item for item in odd_spawn.candidate_evaluations if item.admitted], key=lambda item: (-item.residual_capture, item.candidate.canonical_key()))[0].candidate.canonical_key(),
        "even_selected_candidate_is_highest_score": even_spawn.selected_candidate.canonical_key() == sorted([item for item in even_spawn.candidate_evaluations if item.admitted], key=lambda item: (-item.residual_capture, item.candidate.canonical_key()))[0].candidate.canonical_key(),
        "odd_spawn_projection_loss_bounded": odd_spawn.projection.relative_projection_loss < thresholds["projection_loss"],
        "even_spawn_projection_loss_bounded": even_spawn.projection.relative_projection_loss < thresholds["projection_loss"],
        "odd_spawn_energy_jump_bounded": abs(odd_spawn.projection.energy_jump_hartree) < thresholds["projection_energy_jump_hartree"],
        "even_spawn_energy_jump_bounded": abs(even_spawn.projection.energy_jump_hartree) < thresholds["projection_energy_jump_hartree"],
        "odd_spawn_fidelity_is_unity": 1.0 - odd_spawn.projection.normalized_fidelity < thresholds["projection_loss"],
        "even_spawn_fidelity_is_unity": 1.0 - even_spawn.projection.normalized_fidelity < thresholds["projection_loss"],
        "spawn_stable_id_is_monotone": odd_spawn.added_packet_id == "g000002",
        "spawn_newborn_age_is_zero": odd_spawn.packet_ages_after[-1] == 0,
        "permutation_selects_same_candidate": odd_spawn.selected_candidate.canonical_key() == permutation_spawn.selected_candidate.canonical_key(),
        "permutation_covariance_bounded": permutation_error < thresholds["covariance_error"],
        "constant_gauge_selects_same_candidate": odd_spawn.selected_candidate.canonical_key() == gauge_spawn.selected_candidate.canonical_key(),
        "constant_gauge_covariance_bounded": gauge_error < thresholds["covariance_error"],
        "constant_gauge_score_invariant": abs(odd_spawn.candidate_evaluations[0].residual_capture - gauge_spawn.candidate_evaluations[0].residual_capture) < thresholds["covariance_error"],
        "activation_step_uses_dormant_mask": activation.active_shape_mask.tolist() == [True, True, False],
        "activation_newborn_amplitude_grows": newborn_population > thresholds["minimum_newborn_population"],
        "activation_dormant_shape_is_fixed": dormant_shape_drift < thresholds["dormant_shape_drift"],
        "activation_norm_is_conserved": abs(activation.norm_change) < thresholds["activation_norm_drift"],
        "activation_nonlinear_residual_bounded": activation.nonlinear_residual_norm < thresholds["activation_residual"],
        "prune_event_accepted": prune.event_kind == "prune",
        "prune_removed_low_population_packet": prune.removed_packet_id == "g000001",
        "prune_projection_loss_bounded": prune.projection.relative_projection_loss < thresholds["prune_loss"],
        "prune_energy_jump_bounded": abs(prune.projection.energy_jump_hartree) < thresholds["prune_energy_jump_hartree"],
        "merge_pair_passes_overlap_gate": merge_overlap >= merge_settings.minimum_merge_overlap,
        "merge_event_accepted": merge.event_kind == "merge",
        "merge_removed_exactly_one_packet": merge.before.ngaussian - merge.after.ngaussian == 1,
        "merge_projection_loss_bounded": merge.projection.relative_projection_loss < thresholds["merge_loss"],
        "merge_energy_jump_bounded": abs(merge.projection.energy_jump_hartree) < thresholds["merge_energy_jump_hartree"],
        "controlled_trajectory_contains_one_spawn": controlled.event_counts["spawn"] == 1,
        "controlled_trajectory_contains_activation_noops": controlled.event_counts["none"] == 2,
        "controlled_trajectory_packet_count_grows_once": controlled.maximum_packet_count == state.ngaussian + 1,
        "controlled_trajectory_norm_is_conserved": controlled.maximum_norm_drift < thresholds["activation_norm_drift"],
        "controlled_trajectory_projection_loss_bounded": controlled.maximum_projection_loss < thresholds["projection_loss"],
        "no_event_threshold_produces_no_topology_change": no_event.event_counts == {"none": 1, "spawn": 0, "prune": 0, "merge": 0},
        "no_event_path_reduces_to_v0252": no_event_error < thresholds["no_event_reduction_error"],
        "maximum_count_gate_produces_noop": maximum_count.event_kind == "none" and maximum_count.reason == "maximum packet count reached",
        "duplicate_candidate_is_rejected": duplicate_evaluation.admitted is False,
        "duplicate_candidate_has_zero_novelty": duplicate_evaluation.novelty < 2.0e-12,
        "duplicate_candidate_rank_gate_fires": "rank-deficient-enlarged-basis" in duplicate_evaluation.rejection_reasons,
        "zero_soc_score_toggle_is_identical": zero_score_error < thresholds["zero_soc_error"],
        "zero_soc_coupling_toggle_is_identical": zero_coupling_error < thresholds["zero_soc_error"],
        "controlled_spawning_claim_is_true": V253_CONTROLLED_BASIS_CLAIMS["controlled_residual_driven_spawning_validated"] is True,
        "coefficient_activation_claim_is_true": V253_CONTROLLED_BASIS_CLAIMS["coefficient_only_newborn_activation_validated"] is True,
        "general_aims_claim_remains_false": V253_CONTROLLED_BASIS_CLAIMS["general_aims_branching_validated"] is False,
        "real_pyscf_trajectory_claim_remains_false": V253_CONTROLLED_BASIS_CLAIMS["real_pyscf_soc_trajectory_admitted"] is False,
    }
    if len(checks) != 60:
        raise AssertionError(f"v0.25.3 validation requires 60 checks, found {len(checks)}.")
    audit = ControlledBasisValidationAuditV253(
        checks={name: bool(value) for name, value in checks.items()},
        metrics=metrics,
        thresholds=thresholds,
        passed=bool(all(checks.values())),
    ).validate()
    return ControlledBasisValidationEvidenceV253(
        odd_spawn_event=odd_spawn,
        even_spawn_event=even_spawn,
        permutation_spawn_event=permutation_spawn,
        gauge_spawn_event=gauge_spawn,
        activation_step=activation,
        prune_event=prune,
        merge_event=merge,
        controlled_trajectory=controlled,
        no_event_trajectory=no_event,
        no_event_reference=no_event_reference,
        maximum_count_event=maximum_count,
        duplicate_candidate_evaluation=duplicate_evaluation,
        zero_soc_enabled_evaluation=zero_enabled,
        zero_soc_disabled_evaluation=zero_disabled,
        audit=audit,
    ).validate()


def save_controlled_basis_validation_evidence_v253(path):
    import json
    from pathlib import Path

    evidence = run_controlled_basis_validation_evidence_v253()
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(
        json.dumps(evidence.as_dict(), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return target
