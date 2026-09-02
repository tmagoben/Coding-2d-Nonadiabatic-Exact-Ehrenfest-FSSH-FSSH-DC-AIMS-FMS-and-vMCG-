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

mass = 1.0
omega = 1.0
sigma = 1.0 / np.sqrt(2.0)
q0 = -1.0
p0 = 0.7

x, dx = uniform_grid(-12.0, 12.0, 2048)

alpha = 1.0 / (2.0 * sigma**2)
psi0 = frozen_gaussian(x, q0, p0, alpha)

V = lambda x: harmonic(x, mass=mass, omega=omega)
dV = lambda x: harmonic_gradient(x, mass=mass, omega=omega)
ddV = lambda x: harmonic_hessian(x, mass=mass, omega=omega)

dt = 0.002
steps = 1000

exact = run_split_operator(
    psi0, x, dx, V, mass=mass, dt=dt, steps=steps, store_every=steps
)

tga = run_thawed_gaussian(
    q0, p0, sigma, mass, V, dV, ddV,
    dt=dt, steps=steps, x=x, store_every=steps
)

psi_exact = exact["psi"][-1]
psi_tga = tga["psi"][-1]

fidelity = abs(inner_product(psi_exact, psi_tga, dx)) ** 2

print("Harmonic oscillator")
print("final exact norm =", exact["norm"][-1])
print("final TGA norm   =", (np.vdot(psi_tga, psi_tga) * dx).real)
print("fidelity         =", fidelity)
print("TGA center q     =", tga["q"][-1])
