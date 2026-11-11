from dataclasses import replace
import json

import numpy as np
import pytest

from gaussian_dynamics.analytic_soc_models_v220 import AnalyticSingletTripletSOCProviderV220
from gaussian_dynamics.complex_gauge_v21 import GaugeTransformedOperatorProviderV21, PhaseMixingGaugeV21, random_unitary_v21
from gaussian_dynamics.correlated_gaussian_tdvp_v270 import CorrelatedGaussianSpinorStateV270
from gaussian_dynamics.curved_basis_adaptation_v293 import build_tdvp_snapshot_v293, evaluate_curved_state_v293
from gaussian_dynamics.molecular_frame_bridge_v300 import MolecularFrameBridgeV300, _hash_v300
from gaussian_dynamics.molecular_hamiltonian_v301 import build_molecular_grid_hamiltonian_v301, build_molecular_tdvp_snapshot_v301
from gaussian_dynamics.molecular_midpoint_v302 import (
    MolecularMidpointSettingsV302, _state_payload, _state_from_payload,
    molecular_midpoint_step_v302, run_molecular_midpoint_v302,
    serialize_molecular_trajectory_v302, replay_molecular_trajectory_v302,
)
from gaussian_dynamics.v302_benchmark import free_packet_fixture_v302, run_v302_midpoint_campaign


@pytest.fixture(scope="module")
def free_model():
    return free_packet_fixture_v302()


@pytest.fixture(scope="module")
def one_step(free_model):
    return molecular_midpoint_step_v302(*free_model, 0.04)


def test_temporal_order_reversal_replay_and_conservation():
    evidence = run_v302_midpoint_campaign()
    assert evidence["passed"], evidence
    # A smaller time step cannot remove the finite-grid Gaussian projection defect.
    assert evidence["refinement"][-1]["grid_exact_wavefunction_error"] > 10*evidence["successive_wavefunction_differences"][-1]
    assert evidence["live_molecular_soc_validated"] is False


def test_step_evolves_every_parameter_block_and_preserves_input(free_model, one_step):
    state, _, _ = free_model
    assert one_step.receipt["before"] == _state_payload(state)
    assert one_step.after.time_au == 0.04
    for name in ("q", "p", "width_matrices", "chirp_matrices", "coefficients"):
        assert np.linalg.norm(getattr(one_step.after, name)-getattr(state, name)) > 1e-9
    assert one_step.receipt["scaled_midpoint_residual"] <= 1
    assert all(a["audit"]["passed"] for a in one_step.receipt["center_transport_audits"])


@pytest.mark.parametrize("dt", [0, np.nan, np.inf, -np.inf, True, np.bool_(False)])
def test_invalid_time_steps_fail_closed(free_model, dt):
    with pytest.raises(ValueError, match="dt"):
        molecular_midpoint_step_v302(*free_model, dt)


@pytest.mark.parametrize("field,value", [
    ("absolute_tolerance", 0), ("relative_tolerance", np.nan),
    ("tangent_step", -1), ("tangent_step", np.bool_(True)),
    ("maximum_iterations", True), ("maximum_iterations", 0), ("maximum_iterations", 1.5),
])
def test_invalid_solver_settings_fail_closed(field, value):
    with pytest.raises(ValueError):
        MolecularMidpointSettingsV302(**{field: value}).validate()


@pytest.mark.parametrize("steps", [0, -1, True, 1.5])
def test_invalid_trajectory_lengths_fail_closed(free_model, steps):
    with pytest.raises(ValueError, match="steps"):
        run_molecular_midpoint_v302(*free_model, dt=0.04, steps=steps)


def test_nonconvergence_raises_without_advancing_input(free_model):
    state, bridge, grid = free_model
    before = _state_payload(state)
    with pytest.raises(ValueError, match="failed to converge"):
        molecular_midpoint_step_v302(state, bridge, grid, 0.08,
            settings=MolecularMidpointSettingsV302(maximum_iterations=1))
    assert _state_payload(state) == before


def test_tampered_receipts_and_states_are_rejected(free_model, one_step):
    _, bridge, grid = free_model
    serialized = serialize_molecular_trajectory_v302([one_step])
    payload = json.loads(serialized)
    payload["steps"][0]["dt"] *= 2
    with pytest.raises(ValueError, match="fingerprint"):
        replay_molecular_trajectory_v302(json.dumps(payload), bridge, grid)
    corrupted = replace(one_step, receipt={**one_step.receipt, "scaled_midpoint_residual": 99})
    with pytest.raises(ValueError, match="fingerprint"):
        serialize_molecular_trajectory_v302([corrupted])
    changed = replace(one_step.after, p=one_step.after.p+0.1)
    with pytest.raises(ValueError, match="no longer matches"):
        serialize_molecular_trajectory_v302([replace(one_step, after=changed)])
    with pytest.raises(ValueError, match="contiguous"):
        serialize_molecular_trajectory_v302([one_step, one_step])
    with pytest.raises(ValueError, match="empty"):
        serialize_molecular_trajectory_v302([])


def test_resigned_wrong_result_still_fails_recomputed_replay(free_model, one_step):
    _, bridge, grid = free_model
    payload = json.loads(serialize_molecular_trajectory_v302([one_step]))
    record = payload["steps"][0]
    record["after"]["p"][0][0] += 0.1
    record.pop("fingerprint")
    record["fingerprint"] = _hash_v300(record)
    payload.pop("fingerprint")
    payload["fingerprint"] = _hash_v300(payload)
    with pytest.raises(ValueError, match="does not reproduce"):
        replay_molecular_trajectory_v302(json.dumps(payload), bridge, grid)


def test_replay_binds_operator(free_model, one_step):
    _, bridge, grid = free_model
    with pytest.raises(ValueError, match="Hamiltonian fingerprint"):
        replay_molecular_trajectory_v302(serialize_molecular_trajectory_v302([one_step]),
            bridge, replace(grid, fingerprint="different operator"))


def test_state_serialization_preserves_correlated_widths_exactly():
    state = CorrelatedGaussianSpinorStateV270(
        [[0.1, 0.2]], [[0.3, -0.2]], [[[1.2345, 0.3456], [0.3456, 2.789]]],
        [[[0.02, 0.07], [0.07, -0.03]]], [[0.7+0.1j, -0.2j]], time_au=0.17,
    )
    payload = json.loads(json.dumps(_state_payload(state)))
    restored = _state_from_payload(payload)
    for name in ("q", "p", "width_matrices", "chirp_matrices", "coefficients"):
        np.testing.assert_array_equal(getattr(restored, name), getattr(state, name))
    assert restored.time_au == state.time_au


def test_cached_wilson_sections_match_uncached_soc_tangent_solve():
    gauge = PhaseMixingGaugeV21(random_unitary_v21(4, seed=302),
        np.array([[0.3], [-0.2], [0.1], [0.25]]), np.zeros(4))
    provider = GaugeTransformedOperatorProviderV21(AnalyticSingletTripletSOCProviderV220(), gauge)
    bridge = MolecularFrameBridgeV300(provider, [0.0])
    grid = build_molecular_grid_hamiltonian_v301(bridge, (np.linspace(-0.6, 0.6, 5),))
    state = CorrelatedGaussianSpinorStateV270(
        [[0.13]], [[0.2]], [[[1.1]]], [[[0.08]]], [[0.8, 0.3j, 0.1, -0.2j]],
    ).normalized()
    observed = build_molecular_tdvp_snapshot_v301(state, bridge, grid)
    connection = bridge.connection_field(transport_mode="nac_ode")
    psi = evaluate_curved_state_v293(state, connection, grid.quadrature)
    reference = build_tdvp_snapshot_v293(state, connection, grid.quadrature, grid.apply(psi))
    np.testing.assert_allclose(observed.psi, psi, atol=2e-13, rtol=2e-13)
    np.testing.assert_allclose(observed.snapshot.velocity, reference.velocity, atol=2e-9, rtol=2e-8)
    np.testing.assert_allclose(observed.snapshot.residual_weighted, reference.residual_weighted, atol=2e-10)


def test_midpoint_trajectory_accepts_two_gaussians_and_replays(free_model):
    _, bridge, grid = free_model
    state = CorrelatedGaussianSpinorStateV270(
        [[-0.8], [0.9]], [[0.3], [-0.2]], [[[1.4]], [[1.1]]],
        [[[0.08]], [[-0.03]]], [[0.8], [0.2j]],
    ).normalized()
    step = molecular_midpoint_step_v302(state, bridge, grid, 0.005)
    assert step.receipt["midpoint_metric_rank"] == state.parameter_count
    serialized = serialize_molecular_trajectory_v302([step])
    assert serialize_molecular_trajectory_v302(replay_molecular_trajectory_v302(serialized, bridge, grid)) == serialized


def test_short_soc_trajectory_agrees_in_coordinate_dependent_frames():
    provider = AnalyticSingletTripletSOCProviderV220()
    gauge = PhaseMixingGaugeV21(random_unitary_v21(4, seed=3021),
        np.array([[0.3], [-0.2], [0.1], [0.25]]), np.zeros(4))
    fixed = MolecularFrameBridgeV300(provider, [0.0])
    moving = MolecularFrameBridgeV300(GaugeTransformedOperatorProviderV21(provider, gauge), [0.0])
    axis = np.linspace(-0.8, 0.8, 9)
    h_fixed = build_molecular_grid_hamiltonian_v301(fixed, (axis,))
    h_moving = build_molecular_grid_hamiltonian_v301(moving, (axis,))
    state = CorrelatedGaussianSpinorStateV270(
        [[0.13]], [[0.2]], [[[1.1]]], [[[0.08]]], [[0.8, 0.3j, 0.1, -0.2j]],
    ).normalized()
    rotated = replace(state, coefficients=np.array([gauge.matrix(state.q[0]).conj().T @ state.coefficients[0]]))
    reference = molecular_midpoint_step_v302(state, fixed, h_fixed, 0.04)
    transformed = molecular_midpoint_step_v302(rotated, moving, h_moving, 0.04)
    psi = build_molecular_tdvp_snapshot_v301(reference.after, fixed, h_fixed).psi
    psi_rotated = build_molecular_tdvp_snapshot_v301(transformed.after, moving, h_moving).psi
    expected = np.array([gauge.matrix(q).conj().T @ value for q, value in zip(h_fixed.quadrature.points, psi)])
    # Midpoint in nonlinear parameter charts is not exactly covariant at finite
    # dt. This short-step tolerance covers truncation and tangent differentiation.
    np.testing.assert_allclose(psi_rotated, expected, atol=2e-8, rtol=2e-8)
    assert np.linalg.norm(reference.after.coefficients-state.coefficients) > 1e-6
