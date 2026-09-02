import numpy as np

from gaussian_dynamics.gauge_graph import ElectronicGaugeGraph
from gaussian_dynamics.graph_electronic import GraphElectronicRegistry, rotate_coefficients
from gaussian_dynamics.graph_gaussian import (
    GraphGaussianTBF,
    build_static_graph_gaussian_matrices,
)


def random_unitary(rng,n):
    X=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n))
    Q,_=np.linalg.qr(X)
    return Q

rng=np.random.default_rng(11)
frames={
    "t0":np.eye(2,dtype=complex),
    "c":random_unitary(rng,2),
    "t1":random_unitary(rng,2),
}

graph=ElectronicGaugeGraph(2)
graph.add_overlap("t0","c",frames["t0"].conj().T@frames["c"])
graph.add_overlap("t1","c",frames["t1"].conj().T@frames["c"])
graph.add_overlap("t0","t1",frames["t0"].conj().T@frames["t1"])

registry=GraphElectronicRegistry(graph)
Hglobal=np.array([[0.10,0.02],[0.02,0.34]],complex)
Fglobal=np.array([
    [[0.02,0.01],[0.01,-0.01]],
    [[-0.03,0.0],[0.0,0.04]],
],complex)

for node,R in frames.items():
    registry.add_operator_data(
        node,
        R.conj().T@Hglobal@R,
        np.asarray([R.conj().T@F@R for F in Fglobal]),
    )

A=np.eye(2)
basis=[
    GraphGaussianTBF("t0",np.array([-0.7,0.2]),np.array([0.4,0.1]),A,np.array([1.0,0.0])),
    GraphGaussianTBF("t1",np.array([0.6,-0.1]),np.array([-0.1,0.3]),A,np.array([0.3,0.954])/np.linalg.norm([0.3,0.954])),
]

M=20*np.eye(2)
S,H=build_static_graph_gaussian_matrices(
    basis,registry,M,reference_selector=lambda i,j:"c"
)

# Apply unrelated U(2) gauges at all electronic nodes.
gauges={node:random_unitary(rng,2) for node in graph.nodes}
registry2=registry.gauge_transformed(gauges)
basis2=[]
for tbf in basis:
    basis2.append(GraphGaussianTBF(
        tbf.node,tbf.q,tbf.p,tbf.A,
        rotate_coefficients(tbf.electronic_coefficients,gauges[tbf.node])
    ))

S2,H2=build_static_graph_gaussian_matrices(
    basis2,registry2,M,reference_selector=lambda i,j:"c"
)

print("Gauge-covariant graph-Gaussian matrix elements")
print("-----------------------------------------------")
print("max |S-S'|:",np.max(np.abs(S-S2)))
print("max |H-H'|:",np.max(np.abs(H-H2)))
print("Hermiticity residual H:",np.max(np.abs(H-H.conj().T)))
