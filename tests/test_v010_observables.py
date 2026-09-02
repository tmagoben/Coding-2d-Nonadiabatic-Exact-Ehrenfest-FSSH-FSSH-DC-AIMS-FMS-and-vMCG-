import numpy as np

from gaussian_dynamics.dynamic_gauge_graph import (
    IncrementalElectronicGraph,
    AnalyticCI2DFrameProvider,
)
from gaussian_dynamics.dynamic_graph_aims import DynamicGraphTBF
from gaussian_dynamics.electronic_observables import (
    reduced_electronic_density_graph,
    exact_reduced_electronic_density_diabatic,
    rotate_density_to_frame,
    density_matrix_populations,
    reduced_electronic_density_analytic_ci_diabatic,
)


def test_exact_reduced_density_trace_matches_grid_norm():
    psi=np.zeros((4,5,2),complex)
    psi[...,0]=1.0
    dx=0.2
    dy=0.3
    norm=np.sum(np.abs(psi)**2)*dx*dy

    rho=exact_reduced_electronic_density_diabatic(psi,dx,dy)

    assert np.isclose(np.trace(rho).real,norm)
    assert np.allclose(rho[1],0.0)


def test_density_rotation_preserves_trace_and_spectrum():
    rho=np.array([[0.7,0.1j],[-0.1j,0.3]],complex)
    theta=0.4
    U=np.array([
        [np.cos(theta),-np.sin(theta)],
        [np.sin(theta), np.cos(theta)],
    ],complex)

    rotated=rotate_density_to_frame(rho,U)

    assert np.allclose(np.trace(rotated),np.trace(rho))
    assert np.allclose(
        np.linalg.eigvalsh(rotated),
        np.linalg.eigvalsh(rho),
    )


def test_graph_reduced_density_is_normalized_for_one_tbf():
    provider=AnalyticCI2DFrameProvider()
    manager=IncrementalElectronicGraph(2)
    q=np.array([0.7,0.4])
    node=("n",0)
    manager.add_from_provider(node,q,provider)

    tbf=DynamicGraphTBF(
        uid=0,
        state=1,
        q=q,
        p=np.array([0.2,0.1]),
        A=np.eye(2),
        node=node,
    )

    rho=reduced_electronic_density_graph(
        [1.0+0j],
        [tbf],
        manager.registry,
        node,
        normalize=True,
    )

    assert np.allclose(np.trace(rho),1.0)
    assert np.allclose(density_matrix_populations(rho),[0.0,1.0])


def test_analytic_ci_diabatic_density_is_normalized():
    q=np.array([0.7,0.4])
    tbf=DynamicGraphTBF(
        uid=0,
        state=1,
        q=q,
        p=np.array([0.2,0.1]),
        A=np.eye(2),
        node=("unused",0),
    )
    rho=reduced_electronic_density_analytic_ci_diabatic(
        [1.0+0j],[tbf],normalize=True
    )
    assert np.allclose(np.trace(rho),1.0)
    assert np.isclose(np.linalg.matrix_rank(rho,tol=1e-10),1)
