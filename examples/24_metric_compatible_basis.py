import numpy as np
from gaussian_dynamics.moving_graph_gaussian import metric_compatible_basis_connection,basis_connection_residual
S0=np.array([[1.0,0.2+0.1j],[0.2-0.1j,1.0]],complex)
S1=np.array([[1.0,0.205+0.098j],[0.205-0.098j,1.0]],complex)
seed=np.array([[0.0+0.03j,0.02-0.01j],[-0.01-0.005j,0.0-0.02j]])
dt=0.01
T=metric_compatible_basis_connection(S0,S1,dt,seed)
print('Metric-compatible T:')
print(T)
print('Residual ||Sdot-T-T^dagger||_F =',basis_connection_residual(S0,S1,T,dt))
