import numpy as np

from gaussian_dynamics import (
    CIPassageConfig,
    compare_managed_exact_diabatic_density,
)

out=compare_managed_exact_diabatic_density(
    CIPassageConfig(),
    managed_dt=0.005,
    exact_dt=0.0025,
    exact_grid_n=64,
    spa_order=0,
    spawn_action_threshold=2e-4,
    max_basis=4,
    overlap_block=0.90,
)

print("v0.10 exact versus graph-Gaussian reduced electronic density")
print("-----------------------------------------------------------")
print("exact diabatic populations:  ",out["populations_exact"])
print("managed diabatic populations:",out["populations_managed"])
print("population L2 error:          ",out["population_l2_error"])
print("density Frobenius error:      ",out["density_frobenius_error"])
print("exact purity:                 ",out["purity_exact"])
print("managed purity:               ",out["purity_managed"])
print("exact linear entropy:         ",out["linear_entropy_exact"])
print("managed linear entropy:       ",out["linear_entropy_managed"])

print("\nA norm-conserving Gaussian run may still fail this observable test.")
