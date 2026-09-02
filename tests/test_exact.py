import numpy as np

from gaussian_dynamics import uniform_grid, frozen_gaussian, run_split_operator
from gaussian_dynamics.potentials import harmonic


def test_split_operator_preserves_norm():
    x, dx = uniform_grid(-12.0, 12.0, 2048)

    psi0 = frozen_gaussian(x, q=-1.0, p=0.7, alpha=1.0)
    V = lambda x: harmonic(x, mass=1.0, omega=1.0)

    out = run_split_operator(
        psi0, x, dx, V,
        mass=1.0,
        dt=0.002,
        steps=500,
        store_every=50,
    )

    assert np.max(np.abs(out["norm"] - 1.0)) < 1e-11
