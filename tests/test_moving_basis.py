import numpy as np

from gaussian_dynamics import uniform_grid
from gaussian_dynamics.moving_basis import (
    moving_basis_matrices,
    basis_functions,
    run_moving_gaussian_basis,
)
from gaussian_dynamics.potentials import harmonic, harmonic_gradient


def test_moving_basis_H_is_hermitian():
    x, dx = uniform_grid(-10.0, 10.0, 2048)

    V = lambda x: harmonic(x, mass=1.0, omega=1.0)
    dV = lambda x: harmonic_gradient(x, mass=1.0, omega=1.0)

    S, H, tau, _, _ = moving_basis_matrices(
        x, dx,
        q=np.array([-1.2, 0.4]),
        p=np.array([0.8, -0.3]),
        alpha=1.0,
        mass=1.0,
        potential=V,
        gradient=dV,
    )

    assert np.allclose(S, S.conj().T, atol=1e-12)
    assert np.allclose(H, H.conj().T, atol=1e-11)


def test_overlap_derivative_identity():
    x, dx = uniform_grid(-10.0, 10.0, 2048)

    V = lambda x: harmonic(x, mass=1.0, omega=1.0)
    dV = lambda x: harmonic_gradient(x, mass=1.0, omega=1.0)

    q = np.array([-1.2, 0.4])
    p = np.array([0.8, -0.3])

    S, _, tau, qdot, pdot = moving_basis_matrices(
        x, dx, q, p, 1.0, 1.0, V, dV
    )

    h = 1e-6
    Gp = basis_functions(x, q + h*qdot, p + h*pdot, 1.0)
    Gm = basis_functions(x, q - h*qdot, p - h*pdot, 1.0)

    Sp = Gp.conj().T @ Gp * dx
    Sm = Gm.conj().T @ Gm * dx
    Sdot_fd = (Sp - Sm) / (2*h)

    assert np.allclose(Sdot_fd, tau.conj().T + tau, atol=2e-7)


def test_short_moving_basis_run_preserves_norm():
    x, dx = uniform_grid(-10.0, 10.0, 1024)
    V = lambda x: harmonic(x, mass=1.0, omega=1.0)
    dV = lambda x: harmonic_gradient(x, mass=1.0, omega=1.0)

    out = run_moving_gaussian_basis(
        x, dx,
        q0=[-1.5, -0.3],
        p0=[0.8, 0.2],
        C0=[1.0 + 0j, 0.25 + 0.1j],
        alpha=1.0,
        mass=1.0,
        potential=V,
        gradient=dV,
        dt=0.0005,
        steps=100,
        store_every=10,
    )

    assert np.max(np.abs(out["norm"] - 1.0)) < 2e-8
