import numpy as np
from gaussian_dynamics.moving_graph_gaussian import metric_compatible_basis_connection,basis_connection_residual,moving_basis_coefficient_step,generalized_norm

def test_metric_connection_identity_and_seed_antihermitian_part():
    S0=np.array([[1,.2+.1j],[.2-.1j,1]],complex);S1=np.array([[1,.205+.098j],[.205-.098j,1]],complex);seed=np.array([[.03j,.02-.01j],[-.01-.005j,-.02j]]);dt=.01
    T=metric_compatible_basis_connection(S0,S1,dt,seed);assert basis_connection_residual(S0,S1,T,dt)<1e-13
    anti=lambda A:.5*(A-A.conj().T);assert np.allclose(anti(T),anti(seed),atol=1e-13)

def test_moving_basis_step_preserves_metric_norm():
    S0=np.array([[1,.1],[.1,1]],complex);S1=np.array([[1,.1002],[.1002,1]],complex);H0=np.array([[.1,.01],[.01,.2]],complex);H1=np.array([[.1001,.0101],[.0101,.1999]],complex);dt=1e-3;T=metric_compatible_basis_connection(S0,S1,dt)
    C=np.array([1,.2+.1j],complex);C/=np.sqrt(generalized_norm(C,S0));C1=moving_basis_coefficient_step(C,S0,H0,S1,H1,T,dt);assert abs(generalized_norm(C1,S1)-1)<2e-8
