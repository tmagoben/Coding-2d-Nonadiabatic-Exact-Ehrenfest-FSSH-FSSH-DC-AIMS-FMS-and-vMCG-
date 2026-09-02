import numpy as np

from gaussian_dynamics.gauge_graph import ElectronicGaugeGraph
from gaussian_dynamics.graph_electronic import GraphElectronicRegistry, rotate_coefficients
from gaussian_dynamics.graph_gaussian import (
    GraphGaussianTBF,
    build_static_graph_gaussian_matrices,
    generalized_cayley_step,
    generalized_norm,
)


def random_unitary(rng,n):
    X=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n))
    Q,_=np.linalg.qr(X)
    return Q


def build_problem():
    rng=np.random.default_rng(2)
    frames={"a":np.eye(2,dtype=complex),"b":random_unitary(rng,2),"c":random_unitary(rng,2)}

    graph=ElectronicGaugeGraph(2)
    for u,v in [("a","b"),("b","c"),("a","c")]:
        graph.add_overlap(u,v,frames[u].conj().T@frames[v])

    reg=GraphElectronicRegistry(graph)
    Hglobal=np.array([[0.1,0.025],[0.025,0.35]],dtype=complex)
    Fglobal=np.array([np.diag([0.02,-0.01]),np.diag([-0.03,0.04])],dtype=complex)
    for node,R in frames.items():
        reg.add_operator_data(
            node,R.conj().T@Hglobal@R,
            np.asarray([R.conj().T@F@R for F in Fglobal]),
        )

    A=np.eye(2)
    basis=[
        GraphGaussianTBF("a",np.array([-0.8,0.3]),np.array([0.5,0.1]),A,np.array([1.0,0.0])),
        GraphGaussianTBF("b",np.array([0.0,0.5]),np.array([0.2,-0.1]),A,np.array([0.0,1.0])),
        GraphGaussianTBF("c",np.array([0.7,-0.2]),np.array([-0.1,0.3]),A,np.array([0.7,0.7])/np.sqrt(0.98)),
    ]
    return rng,graph,reg,basis


def test_graph_gaussian_static_matrices_are_hermitian():
    _,_,reg,basis=build_problem()
    M=20*np.eye(2)

    # Use node b as one common local-diabatic reference for every pair.
    S,H=build_static_graph_gaussian_matrices(
        basis,reg,M,reference_selector=lambda i,j:"b"
    )

    assert np.allclose(S,S.conj().T,atol=1e-12)
    assert np.allclose(H,H.conj().T,atol=1e-12)


def test_graph_gaussian_matrices_are_local_gauge_invariant():
    rng,graph,reg,basis=build_problem()
    M=20*np.eye(2)
    S,H=build_static_graph_gaussian_matrices(
        basis,reg,M,reference_selector=lambda i,j:"b"
    )

    gauges={node:random_unitary(rng,2) for node in graph.nodes}
    regp=reg.gauge_transformed(gauges)

    basis_p=[]
    for tbf in basis:
        basis_p.append(GraphGaussianTBF(
            tbf.node,tbf.q,tbf.p,tbf.A,
            rotate_coefficients(tbf.electronic_coefficients,gauges[tbf.node])
        ))

    Sp,Hp=build_static_graph_gaussian_matrices(
        basis_p,regp,M,reference_selector=lambda i,j:"b"
    )

    assert np.allclose(S,Sp,atol=1e-12)
    assert np.allclose(H,Hp,atol=1e-12)


def test_generalized_cayley_preserves_static_graph_basis_norm():
    _,_,reg,basis=build_problem()
    M=20*np.eye(2)
    S,H=build_static_graph_gaussian_matrices(
        basis,reg,M,reference_selector=lambda i,j:"b"
    )

    C=np.array([1.0,0.15+0.05j,-0.1j],complex)
    C/=np.sqrt(generalized_norm(C,S))
    n0=generalized_norm(C,S)

    for _ in range(500):
        C=generalized_cayley_step(C,S,H,dt=0.01)

    n1=generalized_norm(C,S)
    assert abs(n1-n0) < 2e-12
