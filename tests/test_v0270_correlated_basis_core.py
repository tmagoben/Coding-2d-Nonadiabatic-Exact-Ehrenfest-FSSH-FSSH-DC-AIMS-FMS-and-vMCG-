from dataclasses import replace

import numpy as np
import pytest

from gaussian_dynamics.correlated_basis_adaptation_v270 import (
    ControlledCorrelatedBasisSettingsV270,
    CorrelatedSpawnCandidateV270,
    adapt_correlated_basis_once_v270,
    evaluate_correlated_spawn_candidate_v270,
    generate_correlated_spawn_candidates_v270,
    metric_compatible_activation_mask_v270,
    project_correlated_state_v270,
    run_controlled_correlated_dynamics_v270,
)
from gaussian_dynamics.correlated_gaussian_tdvp_v270 import (
    CorrelatedGaussianSpinorStateV270,
    pack_correlated_parameters_v270,
)
from gaussian_dynamics.multidimensional_soc_v260 import two_state_ci_soc_model_v260


def _state():
    return CorrelatedGaussianSpinorStateV270(
        [[-0.25, 0.1]], [[3.0, 0.2]],
        [[[1.7, 0.25], [0.25, 2.6]]],
        [[[0.02, 0.03], [0.03, -0.01]]],
        [[1.0, 0.0]],
    ).normalized()


def test_v0270_candidates_cover_signed_intrinsic_position_and_momentum_axes():
    candidates = generate_correlated_spawn_candidates_v270(_state())
    assert len(candidates) == 8
    assert {item.displacement_kind for item in candidates} == {"position", "momentum"}
    assert {item.coordinate_axis for item in candidates} == {0, 1}
    assert {item.sign for item in candidates} == {-1, 1}
    eigenvalues, eigenvectors = np.linalg.eigh(_state().width_matrices[0])
    for candidate in candidates:
        displacement = (
            candidate.q - _state().q[0]
            if candidate.displacement_kind == "position"
            else candidate.p - _state().p[0]
        )
        direction = eigenvectors[:, candidate.coordinate_axis]
        assert np.linalg.norm(displacement - direction * np.dot(direction, displacement)) < 3.0e-13


def test_v0270_degenerate_principal_axes_fail_closed():
    state = _state()
    isotropic = replace(state, width_matrices=np.asarray([2.0 * np.eye(2)])).normalized()
    near = replace(
        state, width_matrices=np.asarray([np.diag([2.0, 2.0 + 1.0e-10])])
    ).normalized()
    assert generate_correlated_spawn_candidates_v270(isotropic) == ()
    assert generate_correlated_spawn_candidates_v270(near) == ()


def test_v0270_duplicate_candidate_fails_novelty_and_rank_gates():
    state = _state()
    duplicate = CorrelatedSpawnCandidateV270(
        state.q[0], state.p[0], state.width_matrices[0], state.chirp_matrices[0],
        0, "position", 0, 1,
    )
    result = evaluate_correlated_spawn_candidate_v270(
        state, two_state_ci_soc_model_v260(), duplicate
    )
    assert result.admitted is False
    assert result.novelty < 2.0e-12
    assert "rank-deficient-enlarged-basis" in result.rejection_reasons


def test_v0270_spawn_is_exact_and_rotation_covariant():
    state = _state()
    model = two_state_ci_soc_model_v260()
    angle = 0.371
    rotation = np.asarray(
        [[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]]
    )
    event = adapt_correlated_basis_once_v270(state, model)
    transformed = adapt_correlated_basis_once_v270(
        state.coordinate_rotated(rotation), model.coordinate_rotated(rotation)
    )
    assert event.event_kind == transformed.event_kind == "spawn"
    assert event.after.ngaussian == 2
    assert event.added_packet_id == "g000001"
    assert np.max(np.abs(event.after.coefficients[-1])) == 0.0
    assert event.projection.relative_projection_loss == 0.0
    assert np.max(
        np.abs(
            pack_correlated_parameters_v270(event.after.coordinate_rotated(rotation))
            - pack_correlated_parameters_v270(transformed.after)
        )
    ) < 3.0e-8


def test_v0270_correlated_projection_prune_and_merge_pass_gates():
    model = two_state_ci_soc_model_v260(soc_scale=0.0)
    width = np.asarray([[2.0, 0.25], [0.25, 1.4]])
    chirp = np.asarray([[0.02, 0.03], [0.03, -0.01]])
    prune_state = CorrelatedGaussianSpinorStateV270(
        [[-2.0, 0.0], [2.0, 0.0]], [[0.0, 0.0], [0.0, 0.0]],
        [width, width], [chirp, chirp], [[1.0, 0.0], [1.0e-7, 0.0]],
    ).normalized()
    prune = adapt_correlated_basis_once_v270(
        prune_state, model, packet_ids=("g000000", "g000001"),
        packet_ages=(64, 64), next_packet_serial=2,
    )
    merge_state = CorrelatedGaussianSpinorStateV270(
        [[0.0, 0.0], [1.0e-4, 0.0]], [[0.0, 0.0], [0.0, 0.0]],
        [width, width], [chirp, chirp], [[0.7, 0.0], [0.3, 0.0]],
    ).normalized()
    merge = adapt_correlated_basis_once_v270(
        merge_state, model, packet_ids=("g000000", "g000001"),
        packet_ages=(2, 2), next_packet_serial=2,
    )
    assert prune.event_kind == "prune"
    assert prune.projection.relative_projection_loss < 2.0e-7
    assert abs(prune.projection.energy_jump_hartree) < 2.0e-6
    assert merge.event_kind == "merge"
    assert merge.projection.relative_projection_loss < 2.0e-7
    assert abs(merge.projection.energy_jump_hartree) < 2.0e-6


def test_v0270_projection_rejects_vector_or_non_spd_widths():
    state = _state()
    model = two_state_ci_soc_model_v260()
    with pytest.raises(ValueError, match="incompatible"):
        project_correlated_state_v270(
            state, state.q, state.p, [[1.0, 2.0]], [[0.0, 0.0]], model
        )
    with pytest.raises(ValueError, match="positive definite"):
        project_correlated_state_v270(
            state, state.q, state.p, [[[1.0, 0.0], [0.0, -1.0]]],
            state.chirp_matrices, model,
        )


def test_v0270_newborn_activation_controls_full_matrix_block():
    state = _state()
    model = two_state_ci_soc_model_v260()
    dormant = adapt_correlated_basis_once_v270(state, model).after
    assert metric_compatible_activation_mask_v270(
        dormant, model, locked_active_mask=[True, False]
    ).tolist() == [True, False]
    coefficients = dormant.coefficients.copy()
    coefficients[-1, 1] = 0.03
    active = replace(dormant, coefficients=coefficients).normalized()
    assert metric_compatible_activation_mask_v270(
        active, model, locked_active_mask=[True, False]
    ).tolist() == [True, True]


def test_v0270_short_controlled_trajectory_spawns_and_conserves_norm():
    trajectory = run_controlled_correlated_dynamics_v270(
        _state(), two_state_ci_soc_model_v260(), 0.001, 2
    )
    assert trajectory.event_counts == {"none": 1, "spawn": 1, "prune": 0, "merge": 0}
    assert trajectory.maximum_packet_count == 2
    assert trajectory.maximum_norm_drift < 1.0e-8


def test_v0270_packet_cap_and_policy_tampering_fail_closed():
    state = _state()
    model = two_state_ci_soc_model_v260()
    capped = adapt_correlated_basis_once_v270(
        state, model,
        settings=replace(ControlledCorrelatedBasisSettingsV270(), maximum_packet_count=1),
    )
    assert capped.event_kind == "none"
    assert capped.reason == "maximum packet count reached"
    with pytest.raises(ValueError, match="spawn directions"):
        replace(
            ControlledCorrelatedBasisSettingsV270(), spawn_directions="laboratory axes"
        ).validate()
