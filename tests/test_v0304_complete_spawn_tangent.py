"""Independent checks for the recovered complex-coefficient spawning oracle."""
from dataclasses import replace

import numpy as np
import pytest
from scipy.linalg import expm

from gaussian_dynamics.correlated_gaussian_tdvp_v270 import CorrelatedGaussianSpinorStateV270
from gaussian_dynamics.correlated_basis_adaptation_v270 import (
    MultidimensionalSpawnCandidateV270,
    generate_multidimensional_spawn_candidates_v270,
)
from gaussian_dynamics.curved_basis_adaptation_v292 import (
    ConstantConnectionV292,
    CurvedBasisSettingsV292,
    CurvedTDVPSnapshotV292,
    SVDReceiptV292,
    _evaluate_complete_coefficient_tangent_v292,
    actual_linear_residual_reduction_v292,
    candidate_coefficient_block_v292,
    evaluate_curved_spawn_candidate_v292,
    finite_difference_full_jacobian_v292,
    tensor_legendre_quadrature_v292,
)
from gaussian_dynamics.curved_basis_adaptation_v293 import (
    candidate_coefficient_block_v293,
    evaluate_curved_spawn_candidate_v293,
    finite_difference_full_jacobian_v293,
)
from gaussian_dynamics.curved_geometry_v293 import (
    ConnectionFieldV293,
    CoordinateGaugeV293,
    constant_connection_field_v293,
    gauge_transform_connection_v293,
    straight_path_v293,
    dogleg_path_v293,
    wilson_path_v293,
)


SETTINGS = CurvedBasisSettingsV292(
    minimum_predicted_residual_reduction=1e-14,
    minimum_relative_residual_reduction=1e-12,
)


def snapshot(j, residual):
    metric = np.real(j.conj().T @ j)
    s = np.linalg.svd(metric, compute_uv=False)
    kept = s[s > max(1e-12, 1e-10 * s[0])]
    receipt = SVDReceiptV292(s, 1e-12, len(kept), float(kept[0] / kept[-1]), 0, 0)
    return CurvedTDVPSnapshotV292(
        j, metric, np.zeros(j.shape[1]), np.zeros(j.shape[1]), residual, receipt,
    ).validate()


@pytest.fixture(scope="module")
def physical_case():
    state = CorrelatedGaussianSpinorStateV270(
        q=np.array([[-0.25, 0.15]]), p=np.array([[0.12, -0.08]]),
        width_matrices=np.array([[[1.05, 0.14], [0.14, 0.86]]]),
        chirp_matrices=np.array([[[0.05, -0.02], [-0.02, 0.03]]]),
        coefficients=np.array([[0.84 + 0.08j, 0.43 - 0.14j]]),
    ).normalized()
    sy = np.array([[0, -1j], [1j, 0]])
    sz = np.diag([1., -1.])
    generators = np.array([0.19j * sy, 0.13j * sz])
    field = constant_connection_field_v293(generators)
    quad = tensor_legendre_quadrature_v292(((-2., 2.), (-2., 2.)), (4, 4))
    candidate = generate_multidimensional_spawn_candidates_v270(state)[0]
    j = finite_difference_full_jacobian_v293(state, field, quad)
    k = candidate_coefficient_block_v293(candidate, state.nstate, field, quad)
    # A manufactured TDSE defect isolates the newborn imaginary directions.
    # It is a local linear-algebra fixture, not a molecular trajectory.
    raw = 0.13j * k[:, 0] + (0.025 - 0.09j) * k[:, 1]
    jr = np.vstack([j.real, j.imag])
    rr = np.r_[raw.real, raw.imag]
    rr -= jr @ np.linalg.lstsq(jr, rr, rcond=1e-10)[0]
    residual = rr[:len(raw)] + 1j * rr[len(raw):]
    return state, generators, field, quad, candidate, snapshot(j, residual)


def test_pure_imaginary_newborn_direction_is_not_lost(physical_case):
    candidate = physical_case[4]
    j = np.array([[1.], [0.]], dtype=complex)
    k = np.array([[0.], [1.]], dtype=complex)
    snap = snapshot(j, np.array([0., 1j]))
    ev = _evaluate_complete_coefficient_tangent_v292(snap, candidate, k, settings=SETTINGS)
    assert ev.admitted
    assert ev.schur_rank == 2
    assert ev.predicted_residual_reduction == pytest.approx(1.)
    assert actual_linear_residual_reduction_v292(snap, k) == pytest.approx(1.)


@pytest.mark.parametrize("version", [292, 293])
def test_physical_candidate_matches_independent_full_complex_lstsq(physical_case, version):
    state, generators, field, quad, candidate, snap = physical_case
    if version == 292:
        connection = ConstantConnectionV292(generators)
        block = candidate_coefficient_block_v292(candidate, 2, connection, quad)
        j = finite_difference_full_jacobian_v292(state, connection, quad)
        snap = snapshot(j, snap.residual_weighted)
        ev = evaluate_curved_spawn_candidate_v292(state, connection, quad, snap, candidate, settings=SETTINGS)
    else:
        block = candidate_coefficient_block_v293(candidate, 2, field, quad)
        ev = evaluate_curved_spawn_candidate_v293(state, field, quad, snap, candidate, settings=SETTINGS)
    reference = actual_linear_residual_reduction_v292(snap, block)
    assert reference > 1e-5
    assert ev.admitted
    assert ev.predicted_residual_reduction == pytest.approx(reference, rel=3e-6, abs=1e-12)


@pytest.mark.parametrize("coordinate_dependent", [False, True])
def test_recomputed_gauge_change_preserves_spawn_score(physical_case, coordinate_dependent):
    state, _, field, quad, candidate, snap = physical_case
    generator = 0.7j * np.array([[0., 1.], [1., 0.]])
    def angle(r):
        return 0.4 + 0.18 * r[0] - 0.11 * r[1] if coordinate_dependent else 1.
    def matrix(r):
        return expm(angle(r) * generator)
    def derivative(r):
        grad = np.array([0.18, -0.11]) if coordinate_dependent else np.zeros(2)
        return np.array([x * generator @ matrix(r) for x in grad])
    gauge = CoordinateGaugeV293(2, 2, matrix, derivative)
    transformed = gauge_transform_connection_v293(field, gauge)
    changed = replace(state, coefficients=np.array([
        gauge.at(q).conj().T @ c for q, c in zip(state.q, state.coefficients)
    ]))
    # Rebuild the entire Jacobian from the transformed wavefunction; this is
    # not an identity callback or a transformed copy of the score.
    j = finite_difference_full_jacobian_v293(changed, transformed, quad)
    residual = np.array([
        gauge.at(r).conj().T @ value
        for r, value in zip(quad.points, snap.residual_weighted.reshape(-1, 2))
    ]).ravel()
    original = evaluate_curved_spawn_candidate_v293(state, field, quad, snap, candidate, settings=SETTINGS)
    after = evaluate_curved_spawn_candidate_v293(changed, transformed, quad, snapshot(j, residual), candidate, settings=SETTINGS)
    assert after.admitted == original.admitted
    assert after.schur_rank == original.schur_rank
    assert after.predicted_residual_reduction == pytest.approx(original.predicted_residual_reduction, rel=3e-6)
    assert after.relative_residual_reduction == pytest.approx(original.relative_residual_reduction, rel=3e-6)


def test_coordinate_rotation_preserves_recomputed_score(physical_case):
    state, _, field, quad, candidate, snap = physical_case
    theta = 0.63
    rotation = np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])
    changed = replace(state, q=state.q @ rotation.T, p=state.p @ rotation.T,
                      width_matrices=rotation @ state.width_matrices @ rotation.T,
                      chirp_matrices=rotation @ state.chirp_matrices @ rotation.T)
    shifted = replace(candidate, q=rotation @ candidate.q, p=rotation @ candidate.p,
                      width_matrices=rotation @ candidate.width_matrices @ rotation.T,
                      chirp_matrices=rotation @ candidate.chirp_matrices @ rotation.T)
    new_quad = replace(quad, points=quad.points @ rotation.T)
    transformed = ConnectionFieldV293(2, 2, lambda r: np.einsum('am,mij->aij', rotation, field.at(rotation.T @ r)))
    j = finite_difference_full_jacobian_v293(changed, transformed, new_quad)
    original = evaluate_curved_spawn_candidate_v293(state, field, quad, snap, candidate, settings=SETTINGS)
    after = evaluate_curved_spawn_candidate_v293(changed, transformed, new_quad, snapshot(j, snap.residual_weighted), shifted, settings=SETTINGS)
    assert after.admitted == original.admitted
    assert after.predicted_residual_reduction == pytest.approx(original.predicted_residual_reduction, rel=3e-6)


def test_packet_permutation_preserves_recomputed_score(physical_case):
    state, _, field, quad, candidate, _ = physical_case
    two = replace(state, q=np.vstack([state.q, [0.7, -0.6]]),
                  p=np.vstack([state.p, [-0.1, 0.2]]),
                  width_matrices=np.repeat(state.width_matrices, 2, axis=0),
                  chirp_matrices=np.repeat(state.chirp_matrices, 2, axis=0),
                  coefficients=np.vstack([state.coefficients, [0.15j, 0.2]]))
    changed = replace(two, **{name: getattr(two, name)[::-1] for name in
                             ('q', 'p', 'width_matrices', 'chirp_matrices', 'coefficients')})
    j0 = finite_difference_full_jacobian_v293(two, field, quad)
    j1 = finite_difference_full_jacobian_v293(changed, field, quad)
    rng = np.random.default_rng(296)
    residual = rng.normal(size=j0.shape[0]) + 1j * rng.normal(size=j0.shape[0])
    before = evaluate_curved_spawn_candidate_v293(two, field, quad, snapshot(j0, residual), candidate, settings=SETTINGS)
    after = evaluate_curved_spawn_candidate_v293(changed, field, quad, snapshot(j1, residual), replace(candidate, source_packet=1), settings=SETTINGS)
    assert after.admitted == before.admitted
    assert after.predicted_residual_reduction == pytest.approx(before.predicted_residual_reduction, rel=3e-6)


def test_duplicate_complex_coefficient_space_is_rejected(physical_case):
    candidate = physical_case[4]
    k = np.array([[0.], [1.]], dtype=complex)
    j = np.column_stack([k, 1j * k])
    ev = _evaluate_complete_coefficient_tangent_v292(snapshot(j, np.array([1., 0.])), candidate, k, settings=SETTINGS)
    assert not ev.admitted
    assert ev.schur_rank == 0
    assert ev.predicted_residual_reduction == 0.
    assert 'no-new-curved-tangent-direction' in ev.rejection_reasons


@pytest.mark.parametrize("value", [np.nan, np.inf])
def test_nonfinite_defect_is_rejected(physical_case, value):
    snap = physical_case[-1]
    bad = snap.residual_weighted.copy()
    bad[0] = value
    with pytest.raises(ValueError, match='non-finite'):
        replace(snap, residual_weighted=bad).validate()


@pytest.mark.parametrize("path_factory", [straight_path_v293, dogleg_path_v293])
def test_coincident_endpoints_preserve_identity_transport(physical_case, path_factory):
    state, _, field, _, _, _ = physical_case
    point = state.q[0]
    path = path_factory(point, point)
    np.testing.assert_array_equal(path.end, point)
    np.testing.assert_array_equal(wilson_path_v293(field, path), np.eye(2))


def test_quadrature_node_at_packet_center_is_valid(physical_case):
    state, _, field, quad, candidate, _ = physical_case
    points = quad.points.copy()
    points[0] = candidate.q
    block = candidate_coefficient_block_v293(candidate, state.nstate, field, replace(quad, points=points))
    assert np.all(np.isfinite(block))


@pytest.mark.parametrize("matrix", [np.full((2, 2), np.nan), 2 * np.eye(2)])
def test_invalid_segment_transport_fails_closed(physical_case, matrix):
    field = replace(physical_case[2], segment_transport=lambda a, b: matrix)
    with pytest.raises(ValueError, match="segment transport"):
        wilson_path_v293(field, straight_path_v293([0., 0.], [1., 0.]))
