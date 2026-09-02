import numpy as np

from gaussian_dynamics import uniform_grid
from gaussian_dynamics.spawning import (
    TrajectoryBasisFunction,
    should_spawn,
    spawn_child,
    coupled_basis_matrices,
    propagate_static_basis_coefficients,
    state_populations,
)

x, dx = uniform_grid(-10.0, 10.0, 1024)

parent = TrajectoryBasisFunction(
    state=0,
    q=0.0,
    p=0.8,
    alpha=1.0,
)

print("Parent coupling-region test:", should_spawn(parent, threshold=0.005))

child = spawn_child(parent, target_state=1)
basis = [parent, child]

S, H = coupled_basis_matrices(x, dx, basis, mass=1.0)

C = np.array([1.0 + 0.0j, 0.0 + 0.0j])

print("Initial state populations:", state_populations(C, S, basis))

for _ in range(1000):
    C = propagate_static_basis_coefficients(C, S, H, dt=0.05)

print("Final state populations:  ", state_populations(C, S, basis))
print("\nThis example isolates basis growth and coherent interstate coupling.")
print("It is not presented as a full FMS/AIMS trajectory algorithm.")
