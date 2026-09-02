import numpy as np

from gaussian_dynamics import uniform_grid
from gaussian_dynamics.variational import (
    pack_parameters,
    variational_wavefunction,
    tdvp_velocity,
)
from gaussian_dynamics.potentials import harmonic


def test_tdvp_projection_improves_over_zero_velocity():
    x, dx = uniform_grid(-8.0, 8.0, 384)

    theta = pack_parameters(
        coefficients=[1.0 + 0.0j, 0.25 + 0.1j],
        q=[-1.0, 0.3],
        p=[0.8, -0.2],
        alpha=[1.0, 1.2],
        chirp=[0.0, 0.0],
    )

    velocity, info = tdvp_velocity(
        x, dx, theta, mass=1.0, potential=harmonic
    )

    assert np.all(np.isfinite(velocity))
    assert info["residual_norm"] <= info["zero_velocity_residual"] + 1e-12


def test_variational_wavefunction_is_finite():
    x, dx = uniform_grid(-8.0, 8.0, 256)

    theta = pack_parameters(
        coefficients=[1.0 + 0.0j],
        q=[0.0],
        p=[0.5],
        alpha=[1.0],
        chirp=[0.2],
    )

    psi = variational_wavefunction(x, theta)

    assert np.all(np.isfinite(psi))
    assert (np.vdot(psi, psi) * dx).real > 0.0
