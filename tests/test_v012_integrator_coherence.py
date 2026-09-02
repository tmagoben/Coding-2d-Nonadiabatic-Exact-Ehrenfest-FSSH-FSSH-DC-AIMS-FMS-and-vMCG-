import numpy as np

from gaussian_dynamics.moving_basis_v12 import (
    moving_basis_midpoint_cayley_step,
    fixed_basis_cayley_operator,
    endpoint_generalized_norm_error,
)
from gaussian_dynamics.coherence_metrics import (
    offdiagonal_coherence,
    coherence_magnitude,
    coherence_phase_error,
    density_trace_distance,
    bloch_vector,
)


def test_fixed_basis_cayley_is_unitary():
    H=np.array([
        [0.2,0.1+0.03j],
        [0.1-0.03j,-0.1],
    ],complex)
    U=fixed_basis_cayley_operator(H,0.04)

    assert np.allclose(U.conj().T@U,np.eye(2),atol=1e-13)


def test_moving_basis_cayley_reduces_to_fixed_basis_cayley():
    H=np.array([[0.3,0.07],[0.07,-0.2]],complex)
    S=np.eye(2,dtype=complex)
    T=np.zeros((2,2),complex)
    C=np.array([1.0+0j,0.2j])
    C/=np.linalg.norm(C)
    dt=0.05

    a=moving_basis_midpoint_cayley_step(
        C,S,H,S,H,T,dt
    )
    b=fixed_basis_cayley_operator(H,dt)@C

    assert np.allclose(a,b,atol=1e-13)
    assert abs(endpoint_generalized_norm_error(C,S,a,S)) < 1e-13


def test_coherence_metrics_known_density():
    rho=np.array([
        [0.6,0.2*np.exp(1j*0.3)],
        [0.2*np.exp(-1j*0.3),0.4],
    ])

    assert np.isclose(coherence_magnitude(rho),0.2)
    assert np.isclose(abs(offdiagonal_coherence(rho)),0.2)

    rho2=np.array([
        [0.6,0.2*np.exp(1j*0.1)],
        [0.2*np.exp(-1j*0.1),0.4],
    ])

    assert np.isclose(coherence_phase_error(rho,rho2),0.2)
    assert density_trace_distance(rho,rho) < 1e-15

    b=bloch_vector(rho)
    assert b.shape==(3,)
    assert np.isclose(b[2],0.2)
