import numpy as np

from gaussian_dynamics import (
    uniform_grid,
    frozen_gaussian,
    run_split_operator,
    run_thawed_gaussian,
    inner_product,
)
from gaussian_dynamics.potentials import (
    harmonic,
    harmonic_gradient,
    harmonic_hessian,
)


def test_harmonic_tga_matches_classical_center():
    mass = 1.0
    omega = 1.0
    q0 = -1.2
    p0 = 0.6
    sigma = 1.0 / np.sqrt(2.0)

    V = lambda x: harmonic(x, mass=mass, omega=omega)
    dV = lambda x: harmonic_gradient(x, mass=mass, omega=omega)
    ddV = lambda x: harmonic_hessian(x, mass=mass, omega=omega)

    dt = 0.001
    steps = 600
    t = dt * steps

    out = run_thawed_gaussian(
        q0, p0, sigma, mass, V, dV, ddV,
        dt=dt, steps=steps, store_every=steps
    )

    q_exact = q0 * np.cos(omega * t) + p0 / (mass * omega) * np.sin(omega * t)

    assert abs(out["q"][-1] - q_exact) < 1e-10


def test_tga_matches_exact_grid_for_harmonic_oscillator():
    mass = 1.0
    omega = 1.0
    q0 = -1.0
    p0 = 0.7
    sigma = 1.0 / np.sqrt(2.0)

    x, dx = uniform_grid(-12.0, 12.0, 2048)
    alpha = 1.0 / (2.0 * sigma**2)
    psi0 = frozen_gaussian(x, q0, p0, alpha)

    V = lambda x: harmonic(x, mass=mass, omega=omega)
    dV = lambda x: harmonic_gradient(x, mass=mass, omega=omega)
    ddV = lambda x: harmonic_hessian(x, mass=mass, omega=omega)

    dt = 0.002
    steps = 500

    exact = run_split_operator(
        psi0, x, dx, V,
        mass=mass, dt=dt, steps=steps, store_every=steps
    )

    tga = run_thawed_gaussian(
        q0, p0, sigma, mass, V, dV, ddV,
        dt=dt, steps=steps, x=x, store_every=steps
    )

    overlap = inner_product(exact["psi"][-1], tga["psi"][-1], dx)
    fidelity = abs(overlap) ** 2

    assert fidelity > 0.999999
