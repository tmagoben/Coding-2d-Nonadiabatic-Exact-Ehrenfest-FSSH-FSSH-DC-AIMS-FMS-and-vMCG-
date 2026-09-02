import numpy as np

from gaussian_dynamics import (
    uniform_grid,
    frozen_gaussian,
    gaussian_moments,
    analytic_overlap,
    inner_product,
)

x, dx = uniform_grid(-12.0, 12.0, 2048)

alpha = 1.3
g1 = frozen_gaussian(x, q=-1.0, p=0.8, alpha=alpha)
g2 = frozen_gaussian(x, q=0.7, p=-0.3, alpha=alpha)

moments = gaussian_moments(x, g1, dx)

print("Gaussian 1")
print("norm              =", moments["norm"])
print("<x>               =", moments["x_mean"])
print("Var(x)            =", moments["x_variance"])
print("<p>               =", moments["p_mean"])
print("expected Var(x)   =", 1.0 / (2.0 * alpha))

S_numeric = inner_product(g1, g2, dx)
S_analytic = analytic_overlap(-1.0, 0.8, 0.7, -0.3, alpha)

print("\nOverlap")
print("numeric  =", S_numeric)
print("analytic =", S_analytic)
print("error    =", abs(S_numeric - S_analytic))
