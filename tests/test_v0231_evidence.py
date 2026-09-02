from dataclasses import replace

import numpy as np
import pytest

from gaussian_dynamics import (
    ConvergenceLadderObservationV231,
    FrameInvarianceObservationV231,
    IndependentReferenceObservationV231,
    TrackingSpecificationV231,
)


def test_independent_reference_error_is_derived_from_raw_values():
    evidence = IndependentReferenceObservationV231(
        reference_id="reference",
        observable="SOC matrix",
        unit="hartree",
        value_shape=(2,),
        computed_values=(1.0 + 1.0j, 2.0),
        reference_values=(1.0 + 1.0j, 2.0 + 4.0e-6),
        computed_artifact="computed",
        reference_artifact="independent",
        tolerance=1.0e-5,
    )

    assert evidence.error == pytest.approx(4.0e-6)
    assert evidence.passed


def test_reference_rejects_shared_source_and_shape_mismatch():
    reference = IndependentReferenceObservationV231(
        reference_id="reference",
        observable="SOC",
        unit="hartree",
        value_shape=(1,),
        computed_values=(1.0,),
        reference_values=(1.0,),
        computed_artifact="same",
        reference_artifact="same",
        tolerance=1.0e-5,
    )
    with pytest.raises(ValueError, match="must differ"):
        reference.validate()
    with pytest.raises(ValueError, match="size disagrees"):
        replace(
            reference,
            computed_artifact="computed",
            reference_artifact="reference",
            value_shape=(2,),
        ).validate()


def test_convergence_changes_are_derived_between_adjacent_levels():
    ladder = ConvergenceLadderObservationV231(
        kind="basis",
        labels=("small", "medium", "large"),
        observable="SOC norm",
        unit="hartree",
        value_shape=(1,),
        values=((1.0e-3,), (8.0e-4,), (7.95e-4,)),
        source_artifacts=("small-output", "medium-output", "large-output"),
        tolerance=1.0e-5,
    )

    assert ladder.changes == pytest.approx((2.0e-4, 5.0e-6))
    assert ladder.passed


def test_convergence_ladder_rejects_duplicate_sources():
    with pytest.raises(ValueError, match="distinct"):
        ConvergenceLadderObservationV231(
            kind="method",
            labels=("a", "b"),
            observable="SOC norm",
            unit="hartree",
            value_shape=(1,),
            values=((1.0,), (1.0,)),
            source_artifacts=("same", "same"),
            tolerance=1.0e-5,
        ).validate()


def test_frame_residuals_are_derived_and_rotation_must_be_proper():
    frame = FrameInvarianceObservationV231(
        observable="SOC norm",
        unit="hartree",
        value_shape=(1,),
        base_values=(1.0e-3,),
        translated_values=(1.0e-3 + 1.0e-9,),
        rotated_values=(1.0e-3 + 2.0e-9,),
        expected_rotated_values=(1.0e-3,),
        translation_bohr=(0.2, 0.1, -0.3),
        rotation_matrix=((0.0, -1.0, 0.0), (1.0, 0.0, 0.0), (0.0, 0.0, 1.0)),
        source_artifacts=("base", "translated", "rotated"),
        tolerance=1.0e-8,
    )

    assert frame.translation_residual == pytest.approx(1.0e-9)
    assert frame.rotation_residual == pytest.approx(2.0e-9)
    assert frame.passed
    with pytest.raises(ValueError, match="proper rotation"):
        replace(
            frame,
            rotation_matrix=((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, -1.0)),
        ).validate()


def test_tracking_derives_subspace_singular_values_and_margins():
    overlaps = np.tile(np.eye(4, dtype=complex), (3, 3, 1, 1))
    tracking = TrackingSpecificationV231(
        manifold_labels=("first", "second"),
        manifold_state_indices=((0, 1), (2, 3)),
        record_edges=((0, 1), (1, 2)),
        overlap_threshold=0.8,
        margin_threshold=0.1,
    )
    report = tracking.derive(overlaps)

    assert report["minimum_overlap"] == 1.0
    assert report["minimum_margin"] == 1.0
    assert report["passed"]


def test_tracking_rejects_disconnected_graph_and_incomplete_partition():
    tracking = TrackingSpecificationV231(
        manifold_labels=("first", "second"),
        manifold_state_indices=((0, 1), (2,)),
        record_edges=((0, 1), (2, 3)),
        overlap_threshold=0.8,
        margin_threshold=0.1,
    )
    with pytest.raises(ValueError, match="partition"):
        tracking.validate(nrecord=4, nstate=4)
    with pytest.raises(ValueError, match="connected"):
        replace(tracking, manifold_state_indices=((0, 1), (2, 3))).validate(
            nrecord=4, nstate=4
        )


def test_tracking_detects_manifold_leakage():
    overlaps = np.tile(np.eye(4, dtype=complex), (2, 2, 1, 1))
    overlaps[0, 1] = np.asarray(
        [[0.6, 0.0, 0.8, 0.0], [0.0, 0.6, 0.0, 0.8],
         [0.8, 0.0, 0.6, 0.0], [0.0, 0.8, 0.0, 0.6]],
        dtype=complex,
    )
    tracking = TrackingSpecificationV231(
        manifold_labels=("first", "second"),
        manifold_state_indices=((0, 1), (2, 3)),
        record_edges=((0, 1),),
        overlap_threshold=0.8,
        margin_threshold=0.1,
    )

    report = tracking.derive(overlaps)
    assert not report["passed"]
    assert report["minimum_margin"] < 0.0
