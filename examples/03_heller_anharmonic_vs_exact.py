import numpy as np

from gaussian_dynamics import (
    uniform_grid,
    frozen_gaussian,
    run_split_operator,
    run_thawed_gaussian,
    inner_product,
)
from gaussian_dynamics.potentials import (
    quartic,
    quartic_gradient,
    quartic_hessian,
)

mass = 1.0
q0 = -1.5
p0 = 1.2
sigma = 0.55

x, dx = uniform_grid(-12.0, 12.0, 2048)
alpha = 1.0 / (2.0 * sigma**2)
psi0 = frozen_gaussian(x, q0, p0, alpha)

V = quartic
dV = quartic_gradient
ddV = quartic_hessian

dt = 0.001
steps = 1500

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

print("Anharmonic quartic model")
print("fidelity exact vs TGA =", fidelity)
print("exact norm            =", exact["norm"][-1])
print("TGA norm              =", (np.vdot(psi_tga, psi_tga) * dx).real)
print("\nA declining fidelity is expected when the exact packet develops")
print("non-Gaussian structure that one thawed Gaussian cannot represent.")
