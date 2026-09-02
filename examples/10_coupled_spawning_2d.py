import numpy as np

from gaussian_dynamics.spawned_basis_2d import (
    AdiabaticTBF2D,
    midpoint_grid_2d,
    run_coupled_spawned_basis_2d,
)

x,y,X,Y,points,dx,dy=midpoint_grid_2d(-3,3,24,-3,3,24)

parent=AdiabaticTBF2D(
    state=1,
    q=np.array([0.55,0.45]),
    p=np.array([0.6,0.8]),
    A=1.2*np.eye(2),
)

out=run_coupled_spawned_basis_2d(
    points,dx,dy,
    initial_basis=[parent],
    C0=[1.0+0.0j],
    mass=20.0,
    dt=0.0002,
    steps=50,
    spawn_threshold=1e-6,
    overlap_block=0.9,
    max_basis=2,
    store_every=5,
)

print("Coupled spawned Gaussian basis")
print("spawn events:",out["events"])
print("basis-size history:",out["basis_size"])
print("final state populations:",out["state_populations"][-1])
print("max norm drift:",np.max(np.abs(out["norm"]-1.0)))
print("final coefficients:",out["final_coefficients"])
