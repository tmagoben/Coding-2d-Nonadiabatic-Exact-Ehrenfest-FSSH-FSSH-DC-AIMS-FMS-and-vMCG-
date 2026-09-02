import numpy as np

from gaussian_dynamics.dynamic_gauge_graph import (
    IncrementalElectronicGraph,
    AnalyticCI2DFrameProvider,
)
from gaussian_dynamics.graph_gaussian import GraphGaussianTBF
from gaussian_dynamics.spa_matrix_elements_v11 import (
    build_graph_gaussian_matrices_spa_general,
)


def test_unequal_width_graph_spa_matrices_are_hermitian():
    provider=AnalyticCI2DFrameProvider(nuclear_mass_au=20.0)
    manager=IncrementalElectronicGraph(2)

    q0=np.array([-0.9,0.5])
    q1=np.array([0.7,0.8])
    qc=0.5*(q0+q1)

    manager.add_from_provider("a",q0,provider)
    manager.add_from_provider("b",q1,provider,connect_to=["a"])
    manager.add_from_provider("c",qc,provider,connect_to=["a","b"])

    e0=np.array([1.0,0.0],complex)
    e1=np.array([0.0,1.0],complex)

    basis=[
        GraphGaussianTBF(
            "a",
            q0,
            np.array([0.4,0.1]),
            np.array([[1.2,0.1],[0.1,0.8]]),
            e0,
        ),
        GraphGaussianTBF(
            "b",
            q1,
            np.array([-0.2,0.3]),
            np.array([[0.7,-0.05],[-0.05,1.3]]),
            e1,
        ),
    ]

    ref=lambda i,j: "a" if i==j==0 else ("b" if i==j==1 else "c")
    M=20.0*np.eye(2)

    for order in (0,1):
        S,H=build_graph_gaussian_matrices_spa_general(
            basis,
            manager.registry,
            M,
            ref,
            order=order,
        )
        assert np.allclose(S,S.conj().T,atol=2e-12)
        assert np.allclose(H,H.conj().T,atol=2e-11)
