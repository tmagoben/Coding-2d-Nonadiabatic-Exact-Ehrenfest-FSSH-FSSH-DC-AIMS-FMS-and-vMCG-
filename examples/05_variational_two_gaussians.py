import numpy as np

from gaussian_dynamics import (
    uniform_grid,
    pack_parameters,
    tdvp_velocity,
    run_variational_dynamics,
)
from gaussian_dynamics.potentials import quartic

x, dx = uniform_grid(-8.0, 8.0, 384)

theta0 = pack_parameters(
    coefficients=[0.9 + 0.0j, 0.35 + 0.1j],
    q=[-1.3, -0.3],
    p=[1.0, 0.2],
    alpha=[1.2, 0.9],
    chirp=[0.0, 0.0],
)

velocity, info = tdvp_velocity(
    x, dx, theta0, mass=1.0, potential=quartic
)

print("Initial TDVP diagnostics")
print("rank(G)                    =", info["rank"])
print("projected TDSE residual    =", info["residual_norm"])
print("zero-velocity residual     =", info["zero_velocity_residual"])

out = run_variational_dynamics(
    x=x,
    dx=dx,
    theta0=theta0,
    mass=1.0,
    potential=quartic,
    dt=0.0002,
    steps=40,
    store_every=10,
)

print("\nShort variational propagation")
print("initial norm =", out["norm"][0])
print("final norm   =", out["norm"][-1])
print("final residual =", out["residual"][-1])
