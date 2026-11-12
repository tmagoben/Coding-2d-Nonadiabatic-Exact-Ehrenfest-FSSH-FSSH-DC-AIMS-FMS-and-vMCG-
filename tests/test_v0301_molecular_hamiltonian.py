from dataclasses import replace

import numpy as np
import pytest
from scipy.linalg import block_diag, expm

from gaussian_dynamics.analytic_soc_models_v220 import AnalyticSingletTripletSOCProviderV220
from gaussian_dynamics.complex_gauge_v21 import (
    GaugeTransformedOperatorProviderV21, PhaseMixingGaugeV21, random_unitary_v21,
)
from gaussian_dynamics.correlated_gaussian_tdvp_v270 import CorrelatedGaussianSpinorStateV270
from gaussian_dynamics.curved_basis_adaptation_v292 import CurvedQuadratureV292
from gaussian_dynamics.curved_geometry_v293 import (
    constant_connection_field_v293, dogleg_path_v293, straight_path_v293, wilson_path_v293,
)
from gaussian_dynamics.electronic_contract_v213 import (
    ElectronicModelSpaceV213, ElectronicOperatorProvenanceV213, ElectronicStateDescriptorV213,
)
from gaussian_dynamics.molecular_frame_bridge_v300 import (
    MolecularFrameBridgeSettingsV300, MolecularFrameBridgeV300,
)
from gaussian_dynamics.molecular_hamiltonian_v301 import (
    build_molecular_grid_hamiltonian_v301, build_molecular_tdvp_snapshot_v301,
)
from gaussian_dynamics.synthetic_operator_provider_v21 import (
    SyntheticLinearOperatorConfigV21, SyntheticLinearOperatorProviderV21,
)


def scalar_provider(ndim=1, mass=2.0):
    """Explicit zero-potential analytic fixture, not molecular execution evidence."""
    provider = SyntheticLinearOperatorProviderV21(SyntheticLinearOperatorConfigV21(
        nstate=1, nq=ndim, mass=mass, base_scale=0.0, derivative_scale=0.0,
    ))
    space = ElectronicModelSpaceV213(
        "scalar fixture", "fixed_general", (ElectronicStateDescriptorV213("scalar"),),
    )
    provider.provenance = ElectronicOperatorProvenanceV213(
        "zero-potential grid fixture", "1", space, "analytic zero potential",
        parameters={"mass": mass, "ndim": ndim},
    ).validate()
    return provider


def bridge_for(provider, ndim=1):
    return MolecularFrameBridgeV300(
        provider, np.zeros(ndim),
        settings=MolecularFrameBridgeSettingsV300(
            maximum_audit_segment_length=0.6, connection_overlap_tolerance=3e-9,
        ),
    )


@pytest.fixture(scope="module")
def soc_grids():
    provider = AnalyticSingletTripletSOCProviderV220()
    gauge = PhaseMixingGaugeV21(
        random_unitary_v21(4, seed=301),
        np.array([[0.31], [-0.22], [0.17], [0.28]]),
        np.array([0.1, 0.2, -0.1, -0.3]),
    )
    moving = GaugeTransformedOperatorProviderV21(provider, gauge)
    fixed_bridge, moving_bridge = bridge_for(provider), bridge_for(moving)
    axes = (np.linspace(-0.8, 0.8, 9),)
    fixed = build_molecular_grid_hamiltonian_v301(fixed_bridge, axes)
    transformed = build_molecular_grid_hamiltonian_v301(moving_bridge, axes)
    return fixed_bridge, moving_bridge, gauge, fixed, transformed


def test_zero_connection_reduces_to_independent_dirichlet_stencil():
    bridge = bridge_for(scalar_provider())
    axis = np.linspace(-0.8, 0.8, 9)
    grid = build_molecular_grid_hamiltonian_v301(bridge, (axis,))
    h = axis[1] - axis[0]
    expected = (2*np.eye(9) - np.eye(9, k=1) - np.eye(9, k=-1)) / (4*h*h)
    np.testing.assert_allclose(grid.matrix.toarray(), expected, atol=1e-14)
    assert grid.matrix[0, -1] == 0
    assert np.min(np.linalg.eigvalsh(grid.kinetic.toarray())) > 0
    assert len(grid.segment_audits) == 8


def test_two_dimensional_kinetic_matches_independent_kronecker_sum():
    bridge = bridge_for(scalar_provider(ndim=2), ndim=2)
    x, y = np.linspace(-0.4, 0.4, 3), np.linspace(-0.3, 0.3, 4)
    grid = build_molecular_grid_hamiltonian_v301(bridge, (x, y))
    tx = (2*np.eye(3) - np.eye(3, k=1) - np.eye(3, k=-1)) / (4*(x[1]-x[0])**2)
    ty = (2*np.eye(4) - np.eye(4, k=1) - np.eye(4, k=-1)) / (4*(y[1]-y[0])**2)
    np.testing.assert_allclose(grid.matrix.toarray(), np.kron(tx, np.eye(4)) + np.kron(np.eye(3), ty))
    assert len(grid.segment_audits) == 2*4 + 3*3


def test_full_soc_action_is_hermitian_and_pointwise_gauge_covariant(soc_grids):
    _, _, gauge, fixed, moving = soc_grids
    transform = block_diag(*[gauge.matrix(q) for q in fixed.quadrature.points])
    expected = transform.conj().T @ fixed.matrix.toarray() @ transform
    np.testing.assert_allclose(moving.matrix.toarray(), expected, atol=2e-12, rtol=2e-12)
    np.testing.assert_allclose(moving.matrix.toarray(), moving.matrix.toarray().conj().T, atol=2e-14)
    assert np.min(np.linalg.eigvalsh(moving.kinetic.toarray())) > 0
    rng = np.random.default_rng(3011)
    psi = rng.normal(size=(9, 4)) + 1j*rng.normal(size=(9, 4))
    rotated = (transform.conj().T @ psi.reshape(-1)).reshape(psi.shape)
    expected_action = (transform.conj().T @ fixed.apply(psi).reshape(-1)).reshape(psi.shape)
    np.testing.assert_allclose(moving.apply(rotated), expected_action, atol=3e-12)
    electronic_only = np.einsum("nij,nj->ni", fixed.electronic_blocks, psi)
    assert np.linalg.norm(fixed.apply(psi) - electronic_only) > 1e-3


def test_discrete_exact_evolution_preserves_norm_energy_and_gauge(soc_grids):
    _, _, gauge, fixed, moving = soc_grids
    transform = block_diag(*[gauge.matrix(q) for q in fixed.quadrature.points])
    rng = np.random.default_rng(3012)
    psi = rng.normal(size=36) + 1j*rng.normal(size=36)
    psi /= np.linalg.norm(psi)
    evolved = expm(-0.7j * fixed.matrix.toarray()) @ psi
    moved = expm(-0.7j * moving.matrix.toarray()) @ (transform.conj().T @ psi)
    np.testing.assert_allclose(moved, transform.conj().T @ evolved, atol=3e-12)
    assert abs(np.vdot(evolved, evolved).real - 1.0) < 2e-13
    assert abs(np.vdot(evolved, fixed.matrix @ evolved) - np.vdot(psi, fixed.matrix @ psi)) < 2e-13


def test_tdvp_wavefunction_action_and_residual_are_gauge_covariant(soc_grids):
    fixed_bridge, moving_bridge, gauge, fixed, moving = soc_grids
    state = CorrelatedGaussianSpinorStateV270(
        [[0.13]], [[0.2]], [[[1.1]]], [[[0.08]]], [[0.8, 0.3j, 0.1, -0.2j]],
    ).normalized()
    transformed_state = replace(
        state, coefficients=np.array([gauge.matrix(state.q[0]).conj().T @ state.coefficients[0]]),
    )
    reference = build_molecular_tdvp_snapshot_v301(state, fixed_bridge, fixed)
    rotated = build_molecular_tdvp_snapshot_v301(transformed_state, moving_bridge, moving)
    for original, transformed in ((reference.psi, rotated.psi), (reference.hpsi, rotated.hpsi)):
        expected = np.array([gauge.matrix(q).conj().T @ v for q, v in zip(fixed.quadrature.points, original)])
        np.testing.assert_allclose(transformed, expected, atol=3e-11)
    residual = reference.snapshot.residual_weighted.reshape(9, 4)
    expected = np.array([gauge.matrix(q).conj().T @ v for q, v in zip(fixed.quadrature.points, residual)])
    np.testing.assert_allclose(rotated.snapshot.residual_weighted.reshape(9, 4), expected, atol=2e-8)


def test_gaussian_kinetic_action_converges_at_second_order():
    bridge = bridge_for(scalar_provider())
    errors = []
    for spacing in (0.4, 0.2, 0.1):
        axis = np.arange(-8 + spacing, 8 - spacing/2, spacing)
        grid = build_molecular_grid_hamiltonian_v301(bridge, (axis,))
        psi = np.exp(-0.5*axis**2 + 0.3j*axis)
        exact = (1 - (-axis + 0.3j)**2) * psi / 4.0
        observed = grid.apply(psi[:, None])[:, 0]
        errors.append(float(np.linalg.norm(observed - exact) / np.linalg.norm(exact)))
    orders = np.log2(np.array(errors[:-1]) / errors[1:])
    assert np.all(orders > 1.9), (errors, orders)
    assert np.all(orders < 2.1), (errors, orders)


def test_full_hamiltonian_feeds_curved_tdvp_and_recovers_free_gaussian_velocity():
    bridge = bridge_for(scalar_provider())
    grid = build_molecular_grid_hamiltonian_v301(bridge, (np.linspace(-6, 6, 61),))
    # The centre is exactly a grid node: its Wilson line must be identity.
    state = CorrelatedGaussianSpinorStateV270(
        [[0.0]], [[0.3]], [[[1.1]]], [[[0.08]]], [[1.0+0j]],
    ).normalized()
    result = build_molecular_tdvp_snapshot_v301(state, bridge, grid)
    from gaussian_dynamics.correlated_gaussian_tdvp_v270 import _velocity_blocks_v270
    _, qdot, pdot, edot, bdot = _velocity_blocks_v270(state, result.snapshot.velocity)
    np.testing.assert_allclose(qdot, [[0.3/2]], rtol=0.02)
    np.testing.assert_allclose(pdot, [[0.0]], atol=1e-8)
    np.testing.assert_allclose(edot, [[-2*0.08/2]], rtol=0.04)
    np.testing.assert_allclose(bdot, [[(1.1**2 - 0.08**2)/2]], rtol=0.04)
    np.testing.assert_allclose(result.hpsi, grid.apply(result.psi), atol=1e-15)
    assert result.snapshot.solve_receipt.residual_relative < 2e-8
    assert result.hamiltonian_fingerprint == grid.fingerprint


@pytest.mark.parametrize("factory", [straight_path_v293, dogleg_path_v293])
def test_coincident_endpoints_give_identity_without_perturbing_coordinates(factory):
    connection = constant_connection_field_v293(np.array([[[0.4j]]]))
    q = np.array([0.0])
    path = factory(q, q)
    np.testing.assert_array_equal(path.end, q)
    np.testing.assert_array_equal(wilson_path_v293(connection, path), np.eye(1))


@pytest.mark.parametrize("axes", [
    (np.array([0.0]),), (np.array([0.0, np.nan]),), (np.array([0.0, np.inf]),),
    (np.array([0.0, 0.2, 0.5]),), (np.array([0.2, 0.0]),),
    (np.array([0.0, 0.0]),), (np.array([0.0, 0.2j]),),
    (np.array([0.0, 1e-16, 1e-16]),), (np.array([0.0, 1e-16, 3e-16]),),
    (), (np.array([0.0, 0.2]), np.array([0.0, 0.2])),
])
def test_invalid_grids_fail_closed(axes):
    with pytest.raises(ValueError):
        build_molecular_grid_hamiltonian_v301(bridge_for(scalar_provider()), axes)


def test_variable_and_nondiagonal_masses_are_rejected():
    provider = scalar_provider(ndim=2)
    provider.mass_matrix = np.array([[2.0, 0.3], [0.3, 2.0]])
    with pytest.raises(ValueError, match="diagonal mass"):
        build_molecular_grid_hamiltonian_v301(bridge_for(provider, ndim=2), ([-0.1, 0.1], [-0.1, 0.1]))
    base = scalar_provider()

    class VariableMass:
        provenance = base.provenance

        def evaluate_snapshot(self, q):
            snapshot = base.evaluate_snapshot(q)
            return replace(snapshot, point=replace(
                snapshot.point, mass_matrix_q_au=np.array([[2.0 + q[0]**2]]),
            ))

        def snapshot_overlap(self, left, right):
            return base.snapshot_overlap(left, right)

    with pytest.raises(ValueError, match="constant mass"):
        build_molecular_grid_hamiltonian_v301(bridge_for(VariableMass()), ([-0.1, 0.1],))


def test_inconsistent_nac_and_overlap_cannot_enter_kinetic_operator():
    base = scalar_provider()

    class WrongNAC:
        provenance = base.provenance

        def evaluate_snapshot(self, q):
            snapshot = base.evaluate_snapshot(q)
            return replace(snapshot, point=replace(
                snapshot.point, connection_q=np.array([[[0.2j]]]),
            ))

        def snapshot_overlap(self, left, right):
            return base.snapshot_overlap(left, right)

    with pytest.raises(ValueError, match="disagree"):
        build_molecular_grid_hamiltonian_v301(bridge_for(WrongNAC()), ([-0.1, 0.1],))


def test_nonlocal_links_are_rejected():
    with pytest.raises(ValueError, match="local limit"):
        build_molecular_grid_hamiltonian_v301(bridge_for(scalar_provider()), ([-0.4, 0.4],))


def test_action_rejects_nonfinite_values_and_mismatched_quadrature(soc_grids):
    _, _, _, grid, _ = soc_grids
    with pytest.raises(ValueError, match="shape"):
        grid.apply(np.zeros((9, 3)))
    with pytest.raises(ValueError, match="finite"):
        grid.apply(np.full((9, 4), np.nan))
    state = CorrelatedGaussianSpinorStateV270([[0.0]], [[0.0]], [[[1.0]]], [[[0.0]]], [[1, 0, 0, 0]])
    wrong = CurvedQuadratureV292(grid.quadrature.points[::-1], grid.quadrature.weights)
    with pytest.raises(ValueError, match="own grid"):
        grid.hamiltonian_action(state, np.zeros((9, 4)), wrong)


def test_snapshot_rejects_bridge_mismatch_and_bad_difference_step(soc_grids):
    fixed_bridge, moving_bridge, _, grid, _ = soc_grids
    state = CorrelatedGaussianSpinorStateV270([[0.1]], [[0.0]], [[[1.0]]], [[[0.0]]], [[1, 0, 0, 0]])
    with pytest.raises(ValueError, match="fingerprints"):
        build_molecular_tdvp_snapshot_v301(state, moving_bridge, grid)
    for step in (0.0, -1.0, np.nan, np.inf):
        with pytest.raises(ValueError, match="relative_step"):
            build_molecular_tdvp_snapshot_v301(state, fixed_bridge, grid, relative_step=step)


def test_snapshot_rejects_provider_changes_after_grid_assembly():
    provider = scalar_provider()
    bridge = bridge_for(provider)
    grid = build_molecular_grid_hamiltonian_v301(bridge, ([-0.1, 0.1],))
    provider.H0[0, 0] = 0.1
    state = CorrelatedGaussianSpinorStateV270([[0.0]], [[0.0]], [[[1.0]]], [[[0.0]]], [[1.0]])
    with pytest.raises(ValueError, match="changed after"):
        build_molecular_tdvp_snapshot_v301(state, bridge, grid)


def test_operator_fingerprint_replays_and_binds_grid(soc_grids):
    bridge, _, _, grid, _ = soc_grids
    duplicate = build_molecular_grid_hamiltonian_v301(bridge, grid.axes)
    assert duplicate.fingerprint == grid.fingerprint
    changed = build_molecular_grid_hamiltonian_v301(bridge, (np.linspace(-0.8, 0.8, 10),))
    assert changed.fingerprint != grid.fingerprint
    assert grid.evidence()["production_molecular_accuracy_certified"] is False
    assert grid.evidence()["hermiticity_residual"] < 1e-14
