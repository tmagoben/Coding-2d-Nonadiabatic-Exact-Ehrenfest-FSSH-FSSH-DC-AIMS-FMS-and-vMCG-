import numpy as np

from gaussian_dynamics import uniform_grid, run_moving_gaussian_basis
from gaussian_dynamics.potentials import harmonic, harmonic_gradient

x, dx = uniform_grid(-10.0, 10.0, 1024)

mass = 1.0
V = lambda x: harmonic(x, mass=mass, omega=1.0)
dV = lambda x: harmonic_gradient(x, mass=mass, omega=1.0)

out = run_moving_gaussian_basis(
    x=x,
    dx=dx,
    q0=[-1.6, -0.4],
    p0=[1.0, 0.3],
    C0=[1.0 + 0.0j, 0.35 + 0.2j],
    alpha=1.0,
    mass=mass,
    potential=V,
    gradient=dV,
    dt=0.001,
    steps=500,
    store_every=50,
)

print("Moving Gaussian basis")
print("initial norm =", out["norm"][0])
print("final norm   =", out["norm"][-1])
print("max |N-1|   =", np.max(np.abs(out["norm"] - 1.0)))
print("final centers=", out["q"][-1])
print("max cond(S)  =", np.max(out["condition_number"]))
