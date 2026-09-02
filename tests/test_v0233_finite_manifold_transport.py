import numpy as np
import pytest

from gaussian_dynamics.finite_manifold_transport_v233 import (
    FiniteManifoldOverlapPolicyV233,
    analyze_finite_manifold_overlap_v233,
    certified_transport_from_overlap_v233,
    certify_reciprocal_transport_pair_v233,
)
from gaussian_dynamics.overlap_transport import nearest_unitary
from gaussian_dynamics.gauge_graph import (
    ElectronicGaugeGraph,
    nearest_unitary as graph_polar_projection,
)


def _random_unitary(seed, dimension):
    rng = np.random.default_rng(seed)
    raw = rng.normal(size=(dimension, dimension)) + 1j * rng.normal(
        size=(dimension, dimension)
    )
    q, r = np.linalg.qr(raw)
    phases = np.diag(r)
    phases = np.where(np.abs(phases) == 0.0, 1.0, phases / np.abs(phases))
    return q @ np.diag(phases.conj())


def test_physical_contraction_and_unitary_transport_are_distinct():
    overlap = np.diag([0.9, 0.7]).astype(complex)
    result = certified_transport_from_overlap_v233(overlap)

    assert result.physically_consistent
    assert result.trajectory_ready
    assert not np.allclose(result.overlap.conj().T @ result.overlap, np.eye(2))
    assert np.allclose(result.right_to_left_transport, np.eye(2))
    assert np.allclose(
        result.right_to_left_transport.conj().T
        @ result.right_to_left_transport,
        np.eye(2),
    )


def test_polar_transport_is_covariant_under_independent_endpoint_gauges():
    left = _random_unitary(2301, 3)
    right = _random_unitary(2302, 3)
    base = np.diag([0.95, 0.8, 0.65]).astype(complex)
    gauged = left.conj().T @ base @ right

    reference = certified_transport_from_overlap_v233(base)
    transformed = certified_transport_from_overlap_v233(gauged)
    expected = (
        left.conj().T @ reference.right_to_left_transport @ right
    )
    assert np.allclose(transformed.right_to_left_transport, expected, atol=1e-12)
    assert np.allclose(transformed.singular_values, reference.singular_values)


def test_physically_valid_but_unusable_overlap_fails_trajectory_policy():
    policy = FiniteManifoldOverlapPolicyV233(
        minimum_retained_singular_value=0.5,
        maximum_condition_number=10.0,
    )
    overlap = np.diag([0.9, 0.01]).astype(complex)
    report = analyze_finite_manifold_overlap_v233(overlap, policy=policy)

    assert report.physically_consistent
    assert not report.trajectory_ready
    assert "insufficient_manifold_retention" in report.failed_quality_checks
    assert "ill_conditioned_overlap" in report.failed_quality_checks
    with pytest.raises(ValueError, match="not trajectory ready"):
        certified_transport_from_overlap_v233(overlap, policy=policy)


def test_spectral_expansion_and_rank_loss_are_rejected_by_consumers():
    expansion = np.diag([1.01, 0.9]).astype(complex)
    with pytest.raises(ValueError, match="physically inconsistent"):
        certified_transport_from_overlap_v233(expansion)

    rank_lost = np.diag([1.0, 0.0]).astype(complex)
    with pytest.raises(ValueError, match="not trajectory ready"):
        nearest_unitary(rank_lost)


def test_raw_and_transport_reciprocity_are_both_certified():
    left = _random_unitary(2311, 2)
    right = _random_unitary(2312, 2)
    overlap = left @ np.diag([0.9, 0.8]) @ right.conj().T
    pair = certify_reciprocal_transport_pair_v233(
        overlap, overlap.conj().T
    )
    assert pair.overlap_reciprocity_residual < 1e-12
    assert pair.transport_reciprocity_residual < 1e-12

    corrupted = overlap.conj().T.copy()
    corrupted[0, 0] += 1e-3
    with pytest.raises(ValueError, match="adjoint reciprocity"):
        certify_reciprocal_transport_pair_v233(overlap, corrupted)


def test_gauge_optimizer_projection_is_not_misclassified_as_physical_overlap():
    aggregate = np.diag([2.0, 0.75]).astype(complex)
    assert np.allclose(graph_polar_projection(aggregate), np.eye(2))

    graph = ElectronicGaugeGraph(2)
    with pytest.raises(ValueError, match="physically inconsistent"):
        graph.add_overlap("left", "right", aggregate)
