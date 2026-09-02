from dataclasses import replace

import numpy as np
import pytest

from gaussian_dynamics.v230_benchmark import build_v230_reference_replay


def _with_uniform_cross_geometry_contraction(dataset, scale):
    overlaps = np.asarray(dataset.overlaps, dtype=complex).copy()
    nrecord, _, nstate, _ = overlaps.shape
    block = float(scale) * np.eye(nstate, dtype=complex)
    for left in range(nrecord):
        for right in range(left + 1, nrecord):
            overlaps[left, right] = block
            overlaps[right, left] = block.conj().T
    return replace(dataset, overlaps=overlaps)


def test_nonunitary_physical_cross_geometry_overlaps_are_accepted(tmp_path):
    dataset = build_v230_reference_replay(tmp_path)
    contracted = _with_uniform_cross_geometry_contraction(dataset, 0.9)

    assert contracted.validate() is contracted
    diagnostics = contracted.overlap_diagnostics()
    assert diagnostics.maximum_self_identity_residual == 0.0
    assert diagnostics.maximum_reciprocity_residual == 0.0
    assert diagnostics.minimum_cross_geometry_singular_value == pytest.approx(0.9)
    assert diagnostics.maximum_cross_geometry_singular_value == pytest.approx(0.9)
    assert diagnostics.maximum_contraction_excess == 0.0
    assert diagnostics.unordered_cross_geometry_pair_count == (
        len(dataset.q) * (len(dataset.q) - 1) // 2
    )


def test_cross_geometry_spectral_expansion_is_rejected(tmp_path):
    dataset = build_v230_reference_replay(tmp_path)
    overlaps = np.asarray(dataset.overlaps, dtype=complex).copy()
    nstate = overlaps.shape[-1]
    expansive = np.zeros((nstate, nstate), dtype=complex)
    expansive[0, 0] = 0.8
    expansive[0, 1] = 0.8
    assert np.max(np.abs(expansive)) < 1.0
    assert np.linalg.svd(expansive, compute_uv=False)[0] > 1.0
    overlaps[0, 1] = expansive
    overlaps[1, 0] = expansive.conj().T

    corrupted = replace(dataset, overlaps=overlaps)
    with pytest.raises(ValueError, match="expansive; singular values"):
        corrupted.validate()


def test_cross_geometry_reciprocity_corruption_is_rejected(tmp_path):
    dataset = build_v230_reference_replay(tmp_path)
    overlaps = np.asarray(dataset.overlaps, dtype=complex).copy()
    nstate = overlaps.shape[-1]
    overlaps[0, 1] = 0.9 * np.eye(nstate, dtype=complex)
    overlaps[1, 0] = 0.8 * np.eye(nstate, dtype=complex)

    corrupted = replace(dataset, overlaps=overlaps)
    with pytest.raises(ValueError, match="adjoint reciprocity"):
        corrupted.validate()


def test_self_overlap_must_remain_the_identity(tmp_path):
    dataset = build_v230_reference_replay(tmp_path)
    overlaps = np.asarray(dataset.overlaps, dtype=complex).copy()
    overlaps[0, 0] *= 0.999

    corrupted = replace(dataset, overlaps=overlaps)
    with pytest.raises(ValueError, match="exact identity"):
        corrupted.validate()


@pytest.mark.parametrize("tolerance", [-1.0, np.inf, np.nan])
def test_replay_overlap_tolerance_must_be_finite_and_nonnegative(
    tmp_path, tolerance
):
    dataset = build_v230_reference_replay(tmp_path)

    with pytest.raises(ValueError, match="finite and nonnegative"):
        dataset.validate(tolerance=tolerance)
