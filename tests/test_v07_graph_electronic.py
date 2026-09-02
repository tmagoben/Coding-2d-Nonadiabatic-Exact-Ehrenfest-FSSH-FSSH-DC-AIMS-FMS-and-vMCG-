import numpy as np

from gaussian_dynamics.gauge_graph import ElectronicGaugeGraph
from gaussian_dynamics.graph_electronic import (
    derivative_hamiltonian_matrices,
    GraphElectronicRegistry,
    rotate_coefficients,
)


def random_unitary(rng,n):
    X=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n))
    Q,_=np.linalg.qr(X)
    return Q


def test_derivative_hamiltonian_is_hermitian_for_real_adiabatic_data():
    E=np.array([0.1,0.4])
    G=np.array([[1.0,-0.2],[0.3,0.5]])
    D=np.zeros((2,2,2))
    D[0,1]=np.array([0.7,-0.4])
    D[1,0]=-D[0,1]

    F=derivative_hamiltonian_matrices(E,G,D)
    for A in F:
        assert np.allclose(A,A.conj().T,atol=1e-12)


def test_pair_factors_are_invariant_under_arbitrary_local_unitary_gauges():
    rng=np.random.default_rng(22)
    frames={0:np.eye(2,dtype=complex),1:random_unitary(rng,2),2:random_unitary(rng,2)}

    graph=ElectronicGaugeGraph(2)
    graph.add_overlap(0,1,frames[0].conj().T@frames[1])
    graph.add_overlap(1,2,frames[1].conj().T@frames[2])
    graph.add_overlap(0,2,frames[0].conj().T@frames[2])

    reg=GraphElectronicRegistry(graph)

    # Define one physical operator in a global frame, then express it in each
    # node-local frame.
    Hglobal=np.array([[0.2,0.03-0.01j],[0.03+0.01j,0.7]])
    Fglobal=np.array([
        [[0.1,0.02],[0.02,-0.3]],
        [[0.4,-0.01j],[0.01j,0.2]],
    ],dtype=complex)

    for node,R in frames.items():
        reg.add_operator_data(
            node,
            R.conj().T@Hglobal@R,
            np.asarray([R.conj().T@F@R for F in Fglobal]),
        )

    c0=np.array([1.0,0.0],complex)
    c2=np.array([0.2+0.1j,0.9],complex)
    c2/=np.linalg.norm(c2)

    f=reg.pair_factors(0,c0,2,c2,1)

    gauges={node:random_unitary(rng,2) for node in graph.nodes}
    regp=reg.gauge_transformed(gauges)
    c0p=rotate_coefficients(c0,gauges[0])
    c2p=rotate_coefficients(c2,gauges[2])
    fp=regp.pair_factors(0,c0p,2,c2p,1)

    assert np.allclose(f["overlap"],fp["overlap"],atol=1e-12)
    assert np.allclose(f["potential"],fp["potential"],atol=1e-12)
    assert np.allclose(f["derivative_hamiltonian"],
                       fp["derivative_hamiltonian"],atol=1e-12)
