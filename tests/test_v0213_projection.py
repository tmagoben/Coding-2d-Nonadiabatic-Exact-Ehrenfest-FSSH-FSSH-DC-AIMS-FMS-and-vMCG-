import numpy as np
import pytest

from gaussian_dynamics.block_sparse_molecular_v21 import BlockMolecularTBFV21
from gaussian_dynamics.gaussian_nd import gaussian_nd
from gaussian_dynamics.initial_projection_v213 import (
    block_metric_fixed_frame_v213,
    initialize_separable_block_state_v213,
    project_grid_wavefunction_fixed_frame_v213,
    transform_electronic_vector_to_local_frame_v213,
)


def test_four_state_one_dimensional_projection_has_no_two_state_or_2d_assumption():
    x = np.linspace(-9.0, 9.0, 6001)
    points = x[:, None]
    dx = x[1] - x[0]
    basis = [
        BlockMolecularTBFV21(
            0,
            np.asarray([0.35]),
            np.asarray([0.42]),
            np.asarray([[1.4]]),
        )
    ]
    electronic = np.asarray(
        [0.55 + 0.10j, -0.20 + 0.35j, 0.30 - 0.15j, -0.12 - 0.28j]
    )
    electronic /= np.linalg.norm(electronic)
    target = gaussian_nd(points, basis[0].q, basis[0].p, basis[0].A)[..., None] * electronic

    result = project_grid_wavefunction_fixed_frame_v213(
        target, points, dx, basis
    )

    assert result.nstate == 4
    assert result.nuclear_dimension == 1
    assert result.coefficients.shape == (4,)
    assert result.fidelity > 1.0 - 1.0e-13
    assert result.relative_residual < 1.0e-12
    assert np.allclose(result.coefficients, electronic, rtol=0.0, atol=1.0e-11)


def test_projection_and_separable_initializer_support_two_nuclear_dimensions():
    x = np.linspace(-6.0, 6.0, 241)
    y = np.linspace(-5.0, 5.0, 201)
    X, Y = np.meshgrid(x, y, indexing="ij")
    points = np.stack([X, Y], axis=-1)
    weight = (x[1] - x[0]) * (y[1] - y[0])
    basis = [
        BlockMolecularTBFV21(
            3,
            np.asarray([0.2, -0.3]),
            np.asarray([0.1, 0.25]),
            np.asarray([[1.2, 0.15], [0.15, 1.5]]),
        )
    ]
    electronic = np.asarray([1.0, 0.25j, -0.3 + 0.1j], dtype=complex)
    electronic /= np.linalg.norm(electronic)
    target = gaussian_nd(points, basis[0].q, basis[0].p, basis[0].A)[..., None] * electronic

    result = project_grid_wavefunction_fixed_frame_v213(
        target, points, weight, basis
    )
    metric = block_metric_fixed_frame_v213(basis, 3)
    initialized = initialize_separable_block_state_v213(
        np.asarray([2.0 - 0.5j]), electronic, metric
    )

    assert result.nuclear_dimension == 2
    assert result.nstate == 3
    assert result.fidelity > 1.0 - 1.0e-12
    assert result.relative_residual < 1.0e-11
    assert np.isclose(np.real(np.vdot(initialized, metric @ initialized)), 1.0)


def test_explicit_global_to_local_frame_transform():
    frame = np.asarray([[1.0, 1.0j], [1.0j, 1.0]], dtype=complex) / np.sqrt(2.0)
    local = np.asarray([0.2 + 0.4j, -0.1 + 0.3j])
    global_vector = frame @ local
    assert np.allclose(
        transform_electronic_vector_to_local_frame_v213(global_vector, frame),
        local,
        rtol=0.0,
        atol=1.0e-14,
    )
    with pytest.raises(ValueError, match="frame unitarity"):
        transform_electronic_vector_to_local_frame_v213(
            global_vector, 1.01 * frame
        )
