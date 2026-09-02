import numpy as np

from gaussian_dynamics.ci2d import (
    circle_path,
    berry_line_integral,
    parallel_transport_real_state,
    vector_nac_2d,
    adiabatic_energies_2d,
)

path=circle_path(radius=1.2,n=2001)
phase=berry_line_integral(path)
_,final_overlap=parallel_transport_real_state(path,state=0)

print("Two-dimensional conical intersection")
print("------------------------------------")
print("CI energies:", adiabatic_energies_2d(0.0,0.0))
print("d_01 at (1,0.5):",vector_nac_2d((1.0,0.5))[0,1])
print("closed-loop integral of d_01:",phase)
print("pi:",np.pi)
print("transported final/initial eigenvector overlap:",final_overlap)
