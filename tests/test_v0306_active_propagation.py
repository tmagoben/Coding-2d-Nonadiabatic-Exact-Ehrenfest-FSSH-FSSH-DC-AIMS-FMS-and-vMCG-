"""Independent coefficient-space evolution, frozen geometry and activation gates."""
from copy import deepcopy
from dataclasses import replace
import json

import numpy as np
import pytest

from gaussian_dynamics.correlated_gaussian_tdvp_v270 import CorrelatedGaussianSpinorStateV270
from gaussian_dynamics.curved_basis_adaptation_v292 import CurvedBasisSettingsV292
from gaussian_dynamics.curved_basis_adaptation_v293 import active_parameter_indices_v293, coefficient_basis_v293
from gaussian_dynamics.molecular_active_v306 import (
    molecular_active_midpoint_step_v306, activate_molecular_shapes_v306,
    serialize_active_step_v306, replay_active_step_v306,
)
from gaussian_dynamics.molecular_admission_v305 import evaluate_molecular_admission_v305
from gaussian_dynamics.molecular_hamiltonian_v301 import build_molecular_tdvp_snapshot_v301
from gaussian_dynamics.molecular_midpoint_v302 import MolecularMidpointSettingsV302, molecular_midpoint_step_v302, _state_payload
from gaussian_dynamics.molecular_spawn_v304 import _checkpoint_connection
from gaussian_dynamics.pre_ml_contract_v294 import stable_hash_v294
from gaussian_dynamics.v302_benchmark import free_packet_fixture_v302


@pytest.fixture(scope="module")
def model():
    return free_packet_fixture_v302()


@pytest.fixture(scope="module")
def newborn(model):
    return evaluate_molecular_admission_v305(*model, family_id="free", trajectory_id="active", step=0)


@pytest.fixture(scope="module")
def grown(model, newborn):
    return molecular_active_midpoint_step_v306(newborn.after, *model[1:], .02, active_shape_mask=[True, False])


def populated_pair():
    return CorrelatedGaussianSpinorStateV270(
        [[-.8], [.9]], [[.3], [-.2]], [[[1.4]], [[1.1]]],
        [[[.08]], [[-.03]]], [[.8], [.2j]],
    ).normalized()


def test_active_indices_keep_every_complex_coefficient_and_original_order():
    state = populated_pair()
    np.testing.assert_array_equal(active_parameter_indices_v293(state, [True, True]), np.arange(state.parameter_count))
    np.testing.assert_array_equal(active_parameter_indices_v293(state, [True, False]), [0, 1, 2, 3, 4, 6, 8, 10])
    np.testing.assert_array_equal(active_parameter_indices_v293(state, [False, False]), [0, 1, 2, 3])


def test_default_and_all_active_tangents_are_exactly_identical(model):
    a = build_molecular_tdvp_snapshot_v301(*model).snapshot
    b = build_molecular_tdvp_snapshot_v301(*model, active_shape_mask=[True]).snapshot
    for name in ("jacobian", "metric", "rhs", "velocity", "residual_weighted"):
        np.testing.assert_array_equal(getattr(a, name), getattr(b, name))


def test_all_active_midpoint_reduces_exactly_to_v302(model):
    old = molecular_midpoint_step_v302(*model, .005)
    new = molecular_active_midpoint_step_v306(*model, .005, active_shape_mask=[True])
    assert _state_payload(old.after) == _state_payload(new.after)
    assert old.receipt["scaled_midpoint_residual"] == new.receipt["scaled_midpoint_residual"]


def test_frozen_geometry_matches_independent_generalized_crank_nicolson(model):
    state, bridge, grid = model
    state = populated_pair()
    phi = coefficient_basis_v293(state, _checkpoint_connection(bridge), grid.quadrature).reshape(grid.quadrature.npoint, -1)
    gram = phi.conj().T @ (grid.quadrature.weights[:, None] * phi)
    operator = phi.conj().T @ (grid.quadrature.weights[:, None] * (grid.matrix @ phi))
    dt = .01
    expected = np.linalg.solve(gram + .5j*dt*operator, (gram - .5j*dt*operator) @ state.coefficients.ravel())
    step = molecular_active_midpoint_step_v306(state, bridge, grid, dt, active_shape_mask=[False, False])
    np.testing.assert_allclose(step.after.coefficients.ravel(), expected, atol=3e-10, rtol=0)
    for name in ("q", "p", "width_matrices", "chirp_matrices"):
        np.testing.assert_array_equal(getattr(state, name), getattr(step.after, name))
    assert abs(step.receipt["after_observables"]["norm"] - step.receipt["before_observables"]["norm"]) < 2e-10


def test_real_zero_population_newborn_grows_without_moving_shape(newborn, grown):
    assert np.sum(abs(newborn.after.coefficients[-1])**2) < 1e-20
    assert np.sum(abs(grown.after.coefficients[-1])**2) > 1e-6
    for name in ("q", "p", "width_matrices", "chirp_matrices"):
        np.testing.assert_array_equal(getattr(newborn.after, name)[-1], getattr(grown.after, name)[-1])
    assert np.linalg.norm(grown.after.q[0] - newborn.after.q[0]) > 1e-6
    assert grown.receipt["scaled_midpoint_residual"] <= 1
    assert grown.receipt["midpoint_svd_residual"] < 2e-8
    assert grown.receipt["center_transport_audits"][-1]["kind"] == "identity"


def test_active_step_replays_exactly(model, grown):
    serialized = serialize_active_step_v306(grown)
    rebuilt = replay_active_step_v306(serialized, *model[1:])
    assert serialize_active_step_v306(rebuilt) == serialized


@pytest.mark.parametrize("mask", [[1], [], [True, False], ["False"]])
def test_bad_masks_fail_closed(model, mask):
    with pytest.raises(ValueError, match="mask"):
        molecular_active_midpoint_step_v306(*model, .01, active_shape_mask=mask)


@pytest.mark.parametrize("dt", [0, np.nan, np.inf, True])
def test_bad_time_steps_fail_closed(model, dt):
    with pytest.raises(ValueError, match="dt"):
        molecular_active_midpoint_step_v306(*model, dt)


def test_nonconvergence_preserves_source(model):
    before = _state_payload(model[0])
    with pytest.raises(ValueError, match="failed to converge"):
        molecular_active_midpoint_step_v306(*model, .08, settings=MolecularMidpointSettingsV302(maximum_iterations=1))
    assert _state_payload(model[0]) == before


@pytest.mark.parametrize("kind", ["mask", "geometry", "coefficient"])
def test_rehashed_false_evolution_is_rejected(model, grown, kind):
    receipt = deepcopy(grown.receipt)
    if kind == "mask":
        receipt["active_shape_mask"][0] = False
    elif kind == "geometry":
        receipt["after"]["q"][-1][0] += .1
    else:
        receipt["after"]["coefficients_real"][-1][0] += .1
    receipt.pop("fingerprint")
    receipt["fingerprint"] = stable_hash_v294(receipt)
    with pytest.raises(ValueError, match="does not reproduce"):
        replay_active_step_v306(json.dumps(receipt), *model[1:])


def test_zero_coefficient_child_cannot_be_activated(model, newborn):
    mask, receipt = activate_molecular_shapes_v306(newborn.after, *model[1:], [True, False])
    assert mask.tolist() == [True, False]
    assert receipt["trials"][0]["reason"] == "insufficient_coefficient_population"


def test_population_alone_cannot_override_incompatible_metric(model, grown):
    mask, receipt = activate_molecular_shapes_v306(grown.after, *model[1:], [True, False])
    assert receipt["trials"][0]["coefficient_norm_squared"] > 1e-6
    assert not mask[-1]
    assert receipt["trials"][0]["reason"] == "incompatible_trial_metric"


def test_populated_well_conditioned_shape_activates_without_changing_state(model):
    state = populated_pair()
    before = _state_payload(state)
    mask, receipt = activate_molecular_shapes_v306(state, *model[1:], [True, False])
    assert mask.tolist() == [True, True], receipt
    trial = receipt["trials"][0]
    assert trial["projected_shape_rank"] == 4
    assert trial["trial_svd_residual"] < 2e-8
    assert _state_payload(state) == before


@pytest.mark.parametrize("overrides,reason", [
    ({"newborn_activation_population": 1.0}, "insufficient_coefficient_population"),
    ({"maximum_activation_condition_number": 1.0}, "ill_conditioned_new_shape_tangent"),
    ({"maximum_activation_velocity_amplification": 1e-6}, "shape_velocity_amplification"),
])
def test_activation_gates_each_remain_authoritative(model, overrides, reason):
    default = CurvedBasisSettingsV292()
    settings = replace(default, inherited=replace(default.inherited, **overrides))
    mask, receipt = activate_molecular_shapes_v306(populated_pair(), *model[1:], [True, False], settings=settings)
    assert mask.tolist() == [True, False]
    assert receipt["trials"][0]["reason"] == reason
