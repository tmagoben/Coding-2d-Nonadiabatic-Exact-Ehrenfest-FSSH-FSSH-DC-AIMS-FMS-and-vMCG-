from dataclasses import replace
from types import SimpleNamespace

import numpy as np
import pytest

from gaussian_dynamics.adaptive_multigaussian_tdvp_v252 import (
    QuadraticSpinHamiltonianV252,
    ThawedGaussianSpinorStateV252,
    _kinetic_polynomial_v252,
    _tangent_terms_v252,
    adaptive_implicit_midpoint_tdvp_step_v252,
    pack_adaptive_variational_parameters_v252,
    quadratic_spin_hamiltonian_from_provider_v252,
)
from gaussian_dynamics.analytic_soc_models_v220 import (
    AnalyticDoubletSOCProviderV220,
    AnalyticSingletTripletSOCProviderV220,
)
from gaussian_dynamics.complex_gauge_v21 import random_unitary_v21
from gaussian_dynamics.controlled_basis_adaptation_v253 import (
    EVENT_ORDER_V253,
    PROJECTION_POLICY_V253,
    SPAWN_SCORE_V253,
    BasisLifecycleEventV253,
    ControlledBasisSettingsV253,
    CoefficientActivationStepV253,
    SpawnCandidateV253,
    adapt_basis_once_v253,
    build_controlled_metric_system_v253,
    coefficient_activation_implicit_step_v253,
    evaluate_spawn_candidate_v253,
    generate_spawn_candidates_v253,
    project_adaptive_state_v253,
    run_controlled_basis_dynamics_v253,
)


def _state():
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


def _model():
    return quadratic_spin_hamiltonian_from_provider_v252(
        AnalyticDoubletSOCProviderV220()
    )


def _scalar_model():
    return QuadraticSpinHamiltonianV252(
        900.0,
        np.zeros((1, 1)),
        np.zeros((1, 1)),
        np.asarray([[0.004]]),
    ).validate()


def _parameter_error(left, right):
    return float(
        np.max(
            np.abs(
                pack_adaptive_variational_parameters_v252(left)
                - pack_adaptive_variational_parameters_v252(right)
            )
        )
    )


def _packet_values(state, grid):
    result = []
    for q, p, width, chirp in zip(
        state.q, state.p, state.widths, state.chirps
    ):
        y = grid - q
        result.append(
            (width / np.pi) ** 0.25
            * np.exp(-0.5 * width * y**2 + 0.5j * chirp * y**2 + 1.0j * p * y)
        )
    return np.asarray(result)


def test_scope_freezes_residual_projection_order_and_closed_claims():
    settings = ControlledBasisSettingsV253().validate()
    assert "dPsi/dt+iHPsi" in SPAWN_SCORE_V253
    assert "full-SVD" in PROJECTION_POLICY_V253
    assert settings.event_order == EVENT_ORDER_V253
    assert settings.spawning and settings.pruning and settings.merging
    assert settings.shape_activation_population == 1.0e-6
    for name in (
        "multidimensional_nuclear_motion",
        "full_width_matrices",
        "coordinate_dependent_electronic_frame",
        "real_molecular_soc_provider",
        "general_aims_branching",
    ):
        with pytest.raises(ValueError, match="does not admit"):
            replace(settings, **{name: True}).validate()
    with pytest.raises(ValueError, match="one topology event"):
        replace(settings, one_event_per_checkpoint=False).validate()
    with pytest.raises(ValueError, match="spawn, prune, and merge"):
        replace(settings, spawning=False).validate()


@pytest.mark.parametrize(
    "provider",
    [AnalyticDoubletSOCProviderV220(), AnalyticSingletTripletSOCProviderV220()],
)
def test_complete_even_and_odd_spin_models_remain_admitted(provider):
    model = quadratic_spin_hamiltonian_from_provider_v252(provider)
    assert model.nstate == 4
    assert model.complete_spin_manifold is True
    assert model.source["fixed_frame_verified"] is True


def test_candidate_generator_has_four_canonical_directions_per_packet():
    candidates = generate_spawn_candidates_v253(_state())
    assert len(candidates) == 8
    assert {item.displacement_kind for item in candidates} == {
        "position-minus", "position-plus", "momentum-minus", "momentum-plus"
    }
    assert all(item.width > 0.0 for item in candidates)


def test_analytic_candidate_residual_matches_independent_dense_grid():
    state = _state()
    model = _model()
    candidate = generate_spawn_candidates_v253(state)[0]
    receipt = evaluate_spawn_candidate_v253(state, model, candidate)
    grid = np.linspace(-8.0, 8.0, 80001)
    packets = _packet_values(state, grid)
    candidate_packet = _packet_values(
        ThawedGaussianSpinorStateV252(
            q=[candidate.q], p=[candidate.p], widths=[candidate.width],
            chirps=[candidate.chirp], coefficients=[[1.0]],
        ),
        grid,
    )[0]
    dot_wavefunction = np.zeros((state.nstate, len(grid)), dtype=complex)
    for velocity, (packet, vector, polynomial) in zip(
        receipt.metric_velocity, _tangent_terms_v252(state)
    ):
        values = sum(value * grid**degree for degree, value in enumerate(polynomial))
        dot_wavefunction += velocity * vector[:, None] * values[None, :] * packets[packet]
    h_wavefunction = np.zeros_like(dot_wavefunction)
    for packet in range(state.ngaussian):
        kinetic = _kinetic_polynomial_v252(
            state.q[packet], state.p[packet], state.widths[packet],
            state.chirps[packet], model.mass_au,
        )
        kinetic_values = sum(
            value * grid**degree for degree, value in enumerate(kinetic)
        )
        coefficient = state.coefficients[packet]
        for point, x in enumerate(grid):
            potential = model.H0 + x * model.H1 + x * x * model.H2
            h_wavefunction[:, point] += packets[packet, point] * (
                kinetic_values[point] * coefficient + potential @ coefficient
            )
    residual = dot_wavefunction + 1.0j * h_wavefunction
    numerical = np.trapezoid(np.conj(candidate_packet)[None, :] * residual, grid, axis=1)
    assert np.max(np.abs(numerical - receipt.raw_residual_coupling)) < 2.0e-9


def test_duplicate_candidate_is_rejected_by_novelty_and_rank():
    state = _state()
    candidate = SpawnCandidateV253(
        state.q[0], state.p[0], state.widths[0], state.chirps[0], 0, "external"
    )
    receipt = evaluate_spawn_candidate_v253(state, _model(), candidate)
    assert receipt.admitted is False
    assert "insufficient-novelty" in receipt.rejection_reasons
    assert "rank-deficient-enlarged-basis" in receipt.rejection_reasons


def test_spawn_is_exact_projection_with_stable_identity():
    event = adapt_basis_once_v253(
        _state(), _model(),
        packet_ids=("g000000", "g000001"), packet_ages=(2, 2),
        next_packet_serial=2,
    )
    assert event.event_kind == "spawn"
    assert event.added_packet_id == "g000002"
    assert event.packet_ages_after[-1] == 0
    assert event.after.ngaussian == 3
    assert event.projection.relative_projection_loss < 2.0e-13
    assert 1.0 - event.projection.normalized_fidelity < 2.0e-13
    assert abs(event.projection.energy_jump_hartree) < 2.0e-13
    assert np.linalg.norm(event.after.coefficients[-1]) < 2.0e-14


def test_projection_loss_matches_independent_grid_residual():
    state = _state()
    event = adapt_basis_once_v253(
        state, _model(), packet_ids=("a", "b"), packet_ages=(2, 2),
        next_packet_serial=2,
    )
    grid = np.linspace(-8.0, 8.0, 80001)
    source = np.einsum("ig,ia->ag", _packet_values(state, grid), state.coefficients)
    target = np.einsum(
        "ig,ia->ag", _packet_values(event.after, grid), event.after.coefficients
    )
    residual_squared = float(np.real(np.trapezoid(np.sum(np.abs(source - target) ** 2, axis=0), grid)))
    assert residual_squared < 3.0e-12
    assert abs(residual_squared - event.projection.relative_projection_loss) < 3.0e-12


def test_newborn_coefficient_activation_freezes_shape_and_grows_amplitude():
    state = _state()
    model = _model()
    event = adapt_basis_once_v253(
        state, model, packet_ids=("a", "b"), packet_ages=(2, 2),
        next_packet_serial=2,
    )
    step = coefficient_activation_implicit_step_v253(
        event.after, model, 0.02, np.asarray([True, True, False])
    )
    assert isinstance(step, CoefficientActivationStepV253)
    assert step.end.q[-1] == event.after.q[-1]
    assert step.end.p[-1] == event.after.p[-1]
    assert step.end.widths[-1] == event.after.widths[-1]
    assert step.end.chirps[-1] == event.after.chirps[-1]
    assert np.max(np.abs(step.end.q[:-1] - event.after.q[:-1])) > 1.0e-8
    assert np.linalg.norm(step.end.coefficients[-1]) > 1.0e-6
    assert step.nonlinear_residual_norm < 2.0e-10


def test_controlled_reduced_metric_matches_full_metric_on_active_subspace():
    state = _state()
    system = build_controlled_metric_system_v253(
        state, _model(), np.asarray([True, False])
    ).validate()
    indices = system.active_parameter_indices
    assert np.array_equal(
        system.reduced_metric, system.full_metric[np.ix_(indices, indices)]
    )
    assert np.array_equal(system.reduced_rhs, system.full_rhs[indices])
    frozen = np.setdiff1d(np.arange(state.parameter_count), indices)
    assert np.array_equal(system.full_velocity[frozen], np.zeros(len(frozen)))


def test_no_event_path_reduces_to_v0252_implicit_step():
    state = _state()
    model = _model()
    settings = ControlledBasisSettingsV253(
        spawn_residual_capture_threshold=1.0e3
    )
    controlled = run_controlled_basis_dynamics_v253(
        state, model, dt_au=0.02, steps=1, settings=settings
    )
    reference = adaptive_implicit_midpoint_tdvp_step_v252(state, model, 0.02)
    assert controlled.event_counts == {"none": 1, "spawn": 0, "prune": 0, "merge": 0}
    assert _parameter_error(controlled.final_state, reference.end) < 5.0e-15


def test_prune_requires_age_and_accepts_only_small_projection_loss():
    model = _scalar_model()
    state = ThawedGaussianSpinorStateV252(
        q=[-1.0, 1.0], p=[0.0, 0.0], widths=[2.0, 2.0],
        chirps=[0.0, 0.0], coefficients=[[1.0], [1.0e-6]],
    ).normalized()
    young = adapt_basis_once_v253(
        state, model, packet_ids=("a", "b"), packet_ages=(2, 2),
        next_packet_serial=2,
    )
    assert young.event_kind == "none"
    assert "activation" in young.reason
    mature = adapt_basis_once_v253(
        state, model, packet_ids=("a", "b"), packet_ages=(64, 64),
        next_packet_serial=2,
    )
    assert mature.event_kind == "prune"
    assert mature.removed_packet_id == "b"
    assert mature.projection.relative_projection_loss < 2.0e-12
    assert abs(mature.projection.energy_jump_hartree) < 2.0e-8


def test_merge_uses_overlap_then_projection_and_retains_one_geometry():
    model = _scalar_model()
    state = ThawedGaussianSpinorStateV252(
        q=[0.0, 0.01], p=[0.0, 0.0], widths=[2.0, 2.0],
        chirps=[0.0, 0.0], coefficients=[[0.7], [0.3]],
    ).normalized()
    settings = ControlledBasisSettingsV253(
        maximum_merge_projection_loss=1.0e-4,
        maximum_event_energy_jump_hartree=1.0e-4,
    )
    event = adapt_basis_once_v253(
        state, model, packet_ids=("a", "b"), packet_ages=(2, 2),
        next_packet_serial=2, settings=settings,
    )
    assert event.event_kind == "merge"
    assert event.after.ngaussian == 1
    assert event.projection.relative_projection_loss < 1.0e-4
    assert event.after.q[0] in state.q


def test_one_checkpoint_cannot_combine_removal_and_spawn():
    model = _scalar_model()
    state = ThawedGaussianSpinorStateV252(
        q=[-1.0, 1.0], p=[0.0, 0.0], widths=[2.0, 2.0],
        chirps=[0.0, 0.0], coefficients=[[1.0], [1.0e-6]],
    ).normalized()
    event = adapt_basis_once_v253(
        state, model, packet_ids=("a", "b"), packet_ages=(64, 64),
        next_packet_serial=2,
    )
    assert event.event_kind == "prune"
    assert event.removed_packet_id == "b"
    assert event.added_packet_id is None
    assert event.candidate_evaluations == ()


def test_maximum_packet_count_blocks_spawn():
    settings = ControlledBasisSettingsV253(maximum_packet_count=2)
    event = adapt_basis_once_v253(
        _state(), _model(), packet_ids=("a", "b"), packet_ages=(2, 2),
        next_packet_serial=2, settings=settings,
    )
    assert event.event_kind == "none"
    assert event.reason == "maximum packet count reached"


def test_packet_permutation_selects_same_physical_candidate():
    state = _state()
    model = _model()
    direct = adapt_basis_once_v253(
        state, model, packet_ids=("a", "b"), packet_ages=(2, 2),
        next_packet_serial=2,
    )
    order = np.asarray([1, 0])
    permuted = adapt_basis_once_v253(
        state.permuted(order), model, packet_ids=("b", "a"),
        packet_ages=(2, 2), next_packet_serial=2,
    )
    assert direct.selected_candidate.canonical_key() == permuted.selected_candidate.canonical_key()
    restored = permuted.after.permuted(np.asarray([1, 0, 2]))
    assert _parameter_error(direct.after, restored) < 2.0e-12


def test_constant_electronic_gauge_preserves_scores_and_projection():
    state = _state()
    model = _model()
    unitary = random_unitary_v21(4, 25301)
    direct = adapt_basis_once_v253(
        state, model, packet_ids=("a", "b"), packet_ages=(2, 2),
        next_packet_serial=2,
    )
    transformed = adapt_basis_once_v253(
        state.gauge_transformed(unitary), model.gauge_transformed(unitary),
        packet_ids=("a", "b"), packet_ages=(2, 2), next_packet_serial=2,
    )
    assert direct.selected_candidate.canonical_key() == transformed.selected_candidate.canonical_key()
    assert abs(
        direct.candidate_evaluations[0].residual_capture
        - transformed.candidate_evaluations[0].residual_capture
    ) < 2.0e-14
    expected = direct.after.gauge_transformed(unitary)
    assert _parameter_error(transformed.after, expected) < 2.0e-12


def test_short_controlled_trajectory_records_spawn_activation_and_chain():
    trajectory = run_controlled_basis_dynamics_v253(
        _state(), _model(), dt_au=0.02, steps=3
    )
    assert trajectory.event_counts == {"none": 2, "spawn": 1, "prune": 0, "merge": 0}
    assert trajectory.final_state.ngaussian == 3
    assert trajectory.final_packet_ids == ("g000000", "g000001", "g000002")
    assert trajectory.final_packet_ages == (3, 3, 2)
    assert trajectory.maximum_norm_drift < 3.0e-8
    assert trajectory.maximum_projection_loss < 2.0e-10


def test_static_provider_and_rank_deficient_projection_fail_closed():
    static = SimpleNamespace(evaluate_snapshot=lambda q: SimpleNamespace(matrices=1))
    with pytest.raises(TypeError, match="explicit operator provenance"):
        quadratic_spin_hamiltonian_from_provider_v252(static)
    state = _state()
    with pytest.raises(ValueError, match="rank deficient"):
        project_adaptive_state_v253(
            state, _model(), q=[0.0, 0.0], p=[0.0, 0.0],
            widths=[2.0, 2.0], chirps=[0.0, 0.0], event_kind="prune",
        )


def test_projection_and_event_tampering_are_rejected():
    event = adapt_basis_once_v253(
        _state(), _model(), packet_ids=("a", "b"), packet_ages=(2, 2),
        next_packet_serial=2,
    )
    with pytest.raises(ValueError, match="relative_projection_loss"):
        replace(
            event.projection,
            relative_projection_loss=event.projection.relative_projection_loss + 1.0e-3,
        ).validate()
    with pytest.raises(ValueError, match="packet ID"):
        replace(event, packet_ids_after=("a", "b", "wrong")).validate()
    with pytest.raises(ValueError, match="residual_capture"):
        replace(
            event.candidate_evaluations[0],
            residual_capture=event.candidate_evaluations[0].residual_capture + 1.0e-3,
        ).validate()


def test_activation_receipt_tampering_is_rejected():
    event = adapt_basis_once_v253(
        _state(), _model(), packet_ids=("a", "b"), packet_ages=(2, 2),
        next_packet_serial=2,
    )
    step = coefficient_activation_implicit_step_v253(
        event.after, _model(), 0.02, np.asarray([True, True, False])
    )
    with pytest.raises(ValueError, match="inactive newborn shape"):
        replace(step, end=replace(step.end, q=step.end.q + [0.0, 0.0, 1.0e-3])).validate()
    with pytest.raises(ValueError, match="nonlinear residual"):
        replace(step, nonlinear_residual=step.nonlinear_residual + 1.0e-3).validate()
