from dataclasses import replace

import numpy as np
import pytest

from gaussian_dynamics.multidimensional_basis_adaptation_v260 import (
    ControlledMultidimensionalBasisSettingsV260,
    MultidimensionalSpawnCandidateV260,
    adapt_multidimensional_basis_once_v260,
    evaluate_multidimensional_spawn_candidate_v260,
    generate_multidimensional_spawn_candidates_v260,
    metric_compatible_activation_mask_v260,
    run_controlled_multidimensional_dynamics_v260,
)
from gaussian_dynamics.multidimensional_gaussian_tdvp_v260 import (
    DiagonalGaussianSpinorStateV260,
)
from gaussian_dynamics.multidimensional_soc_v260 import two_state_ci_soc_model_v260


def _state():
    return DiagonalGaussianSpinorStateV260(
        q=[[-0.25, 0.0]], p=[[3.0, 0.0]], widths=[[2.0, 2.0]],
        chirps=[[0.0, 0.0]], coefficients=[[1.0, 0.0]],
    ).normalized()


def test_v0260_candidate_set_covers_signed_q_and_p_axes():
    candidates = generate_multidimensional_spawn_candidates_v260(_state())
    assert len(candidates) == 8
    assert {item.displacement_kind for item in candidates} == {"position", "momentum"}
    assert {item.coordinate_axis for item in candidates} == {0, 1}
    assert {item.sign for item in candidates} == {-1, 1}


def test_v0260_duplicate_candidate_fails_novelty_and_rank_gates():
    state = _state()
    candidate = MultidimensionalSpawnCandidateV260(
        state.q[0], state.p[0], state.widths[0], state.chirps[0], 0, "position", 0, 1
    )
    result = evaluate_multidimensional_spawn_candidate_v260(
        state, two_state_ci_soc_model_v260(), candidate
    )
    assert result.admitted is False
    assert result.novelty < 2.0e-12
    assert "rank-deficient-enlarged-basis" in result.rejection_reasons


def test_v0260_spawn_is_exact_and_has_stable_identity():
    event = adapt_multidimensional_basis_once_v260(_state(), two_state_ci_soc_model_v260())
    assert event.event_kind == "spawn"
    assert event.after.ngaussian == 2
    assert event.added_packet_id == "g000001"
    assert event.packet_ages_after[-1] == 0
    assert np.max(np.abs(event.after.coefficients[-1])) == 0.0
    assert event.projection.relative_projection_loss == 0.0


def test_v0260_prune_and_merge_require_projection_gates():
    model = two_state_ci_soc_model_v260(soc_scale=0.0)
    prune_state = DiagonalGaussianSpinorStateV260(
        [[-2.0, 0.0], [2.0, 0.0]], [[0.0, 0.0], [0.0, 0.0]],
        [[2.0, 2.0], [2.0, 2.0]], [[0.0, 0.0], [0.0, 0.0]],
        [[1.0, 0.0], [1.0e-7, 0.0]],
    ).normalized()
    prune = adapt_multidimensional_basis_once_v260(
        prune_state, model, packet_ids=("g000000", "g000001"),
        packet_ages=(64, 64), next_packet_serial=2,
    )
    merge_state = DiagonalGaussianSpinorStateV260(
        [[0.0, 0.0], [1.0e-4, 0.0]], [[0.0, 0.0], [0.0, 0.0]],
        [[2.0, 2.0], [2.0, 2.0]], [[0.0, 0.0], [0.0, 0.0]],
        [[0.7, 0.0], [0.3, 0.0]],
    ).normalized()
    merge = adapt_multidimensional_basis_once_v260(
        merge_state, model, packet_ids=("g000000", "g000001"),
        packet_ages=(2, 2), next_packet_serial=2,
    )
    assert prune.event_kind == "prune"
    assert prune.projection.relative_projection_loss < 2.0e-7
    assert merge.event_kind == "merge"
    assert merge.projection.relative_projection_loss < 2.0e-7


def test_v0260_activation_requires_population_and_metric_safety():
    model = two_state_ci_soc_model_v260()
    spawned = adapt_multidimensional_basis_once_v260(_state(), model).after
    assert metric_compatible_activation_mask_v260(
        spawned, model, locked_active_mask=[True, False]
    ).tolist() == [True, False]
    coefficients = spawned.coefficients.copy()
    # This population is large enough to pass both the frozen condition-number
    # and velocity-amplification gates; 0.01 deliberately remains just outside
    # the latter for this geometry.
    coefficients[-1, 1] = 0.03
    activated = replace(spawned, coefficients=coefficients).normalized()
    assert metric_compatible_activation_mask_v260(
        activated, model, locked_active_mask=[True, False]
    ).tolist() == [True, True]


def test_v0260_controlled_trajectory_spawns_once_and_conserves_norm():
    trajectory = run_controlled_multidimensional_dynamics_v260(
        _state(), two_state_ci_soc_model_v260(), 0.01, 6
    )
    assert trajectory.event_counts == {"none": 5, "spawn": 1, "prune": 0, "merge": 0}
    assert trajectory.maximum_norm_drift < 1.0e-8
    assert trajectory.maximum_packet_count == 2


def test_v0260_packet_cap_and_serial_collision_fail_closed():
    state = _state()
    model = two_state_ci_soc_model_v260()
    capped = adapt_multidimensional_basis_once_v260(
        state, model,
        settings=replace(ControlledMultidimensionalBasisSettingsV260(), maximum_packet_count=1),
    )
    assert capped.event_kind == "none"
    assert capped.reason == "maximum packet count reached"
    with pytest.raises(ValueError, match="collides"):
        adapt_multidimensional_basis_once_v260(
            state, model, packet_ids=("g000001",), packet_ages=(0,), next_packet_serial=1
        )
