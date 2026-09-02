import numpy as np

from gaussian_dynamics import uniform_grid
from gaussian_dynamics.spawning import (
    TrajectoryBasisFunction,
    spawn_child,
    should_spawn,
    coupled_basis_matrices,
    propagate_static_basis_coefficients,
)


def test_spawn_child_changes_state_but_preserves_phase_space_center():
    parent = TrajectoryBasisFunction(state=0, q=0.1, p=0.7, alpha=1.2)
    child = spawn_child(parent, 1)

    assert child.state == 1
    assert child.q == parent.q
    assert child.p == parent.p
    assert child.alpha == parent.alpha


def test_coupling_region_triggers_simple_spawn_rule():
    parent = TrajectoryBasisFunction(state=0, q=0.0, p=0.7, alpha=1.0)
    assert should_spawn(parent, threshold=0.005)


def test_coupled_basis_matrices_are_hermitian_and_cayley_preserves_norm():
    x, dx = uniform_grid(-10.0, 10.0, 2048)

    parent = TrajectoryBasisFunction(state=0, q=0.0, p=0.7, alpha=1.0)
    child = spawn_child(parent, 1)
    basis = [parent, child]

    S, H = coupled_basis_matrices(x, dx, basis, mass=1.0)

    assert np.allclose(S, S.conj().T, atol=1e-12)
    assert np.allclose(H, H.conj().T, atol=1e-11)

    C = np.array([1.0 + 0.0j, 0.0 + 0.0j])
    n0 = np.real(np.vdot(C, S @ C))

    for _ in range(100):
        C = propagate_static_basis_coefficients(C, S, H, dt=0.02)

    n1 = np.real(np.vdot(C, S @ C))

    assert abs(n1 - n0) < 1e-11
