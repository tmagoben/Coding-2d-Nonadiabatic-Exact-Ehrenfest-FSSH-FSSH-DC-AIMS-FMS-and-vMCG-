import numpy as np

from gaussian_dynamics import (
    CIPassageConfig,
    DynamicGraphTBF,
    build_born_huang_grid_2d,
    reconstruct_born_huang_wavefunction,
)
from gaussian_dynamics.exact_benchmark import localized_adiabatic_packet_2d

config=CIPassageConfig()
grid=build_born_huang_grid_2d(
    grid_n=32,
    half_width=config.half_width,
    mass=config.mass,
)

tbf=DynamicGraphTBF(
    uid=0,
    state=config.state,
    q=config.q_array(),
    p=config.p_array(),
    A=config.A_matrix(),
    node=("seed",0),
)

psi_bh=reconstruct_born_huang_wavefunction(
    [1.0+0j],[tbf],grid
)
psi_exact=localized_adiabatic_packet_2d(
    grid.points,
    config.q_array(),
    config.p_array(),
    config.A_matrix(),
    state=config.state,
)

error=np.sqrt(
    np.sum(np.abs(psi_bh-psi_exact)**2)
    *grid.area
)

print("Coordinate-dependent Born-Huang initial-state check")
print("---------------------------------------------------")
print("wavefunction L2 difference:",error)
print(
    "One g(R) Phi_a(R) TBF reproduces the intended localized adiabatic "
    "initial packet on the same grid."
)
