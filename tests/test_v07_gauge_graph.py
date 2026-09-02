import numpy as np

from gaussian_dynamics.gauge_graph import ElectronicGaugeGraph, nearest_unitary
from gaussian_dynamics.ci2d import analytic_adiabatic_vectors


def random_unitary(rng, n):
    X = rng.normal(size=(n,n)) + 1j*rng.normal(size=(n,n))
    Q, R = np.linalg.qr(X)
    phase = np.diag(R)
    phase = np.where(np.abs(phase) > 0, phase/np.abs(phase), 1.0)
    return Q @ np.diag(np.conj(phase))


def test_flat_u2_graph_can_be_globally_trivialized():
    rng = np.random.default_rng(4)
    frames = {0: np.eye(2,dtype=complex)}
    for node in (1,2,3):
        frames[node] = random_unitary(rng,2)

    graph = ElectronicGaugeGraph(2)
    for u,v in [(0,1),(1,2),(2,3),(3,0),(0,2)]:
        graph.add_overlap(u,v,frames[u].conj().T @ frames[v])

    gauges,_ = graph.spanning_tree_gauges(0)
    assert graph.gauge_objective(gauges) < 1e-24

    for cycle in graph.fundamental_cycles(0):
        W = graph.wilson_loop(cycle)
        assert np.allclose(W,np.eye(2),atol=1e-12)


def test_wilson_loop_is_invariant_under_local_gauge_transform():
    rng=np.random.default_rng(8)
    graph=ElectronicGaugeGraph(2)

    U01=random_unitary(rng,2)
    U12=random_unitary(rng,2)
    U20=random_unitary(rng,2)

    graph.add_overlap(0,1,U01)
    graph.add_overlap(1,2,U12)
    graph.add_overlap(2,0,U20)

    W=graph.wilson_loop([0,1,2])

    G={i:random_unitary(rng,2) for i in (0,1,2)}
    transformed=ElectronicGaugeGraph(2)
    for u,v in [(0,1),(1,2),(2,0)]:
        transformed.add_overlap(
            u,v,G[u].conj().T @ graph.overlap(u,v) @ G[v]
        )

    Wp=transformed.wilson_loop([0,1,2])

    assert np.allclose(np.sort_complex(np.linalg.eigvals(W)),
                       np.sort_complex(np.linalg.eigvals(Wp)),atol=1e-12)
    assert np.allclose(np.trace(W),np.trace(Wp),atol=1e-12)


def test_u1_pi_holonomy_cannot_be_gauged_away():
    graph=ElectronicGaugeGraph(1)
    graph.add_overlap(0,1,np.array([[1.0]]))
    graph.add_overlap(1,2,np.array([[1.0]]))
    graph.add_overlap(2,0,np.array([[-1.0]]))

    W=graph.wilson_loop([0,1,2])
    assert np.allclose(W,[[-1.0]])

    tree_gauges,_=graph.spanning_tree_gauges(0)
    tree_obj=graph.gauge_objective(tree_gauges)
    assert tree_obj > 1.0

    sync=graph.synchronize(0,max_iter=100,tolerance=1e-10,restarts=2)
    sync_obj=graph.gauge_objective(sync)

    # Synchronization may distribute the unavoidable phase frustration, but it
    # cannot eliminate the pi Wilson loop.
    assert sync_obj <= tree_obj + 1e-10
    assert np.allclose(graph.wilson_loop([0,1,2]),[[-1.0]])


def test_lower_ci_state_ring_has_minus_one_discrete_holonomy():
    n=120
    phi=np.linspace(-np.pi,np.pi,n,endpoint=False)
    states=[]
    for angle in phi:
        R=np.array([np.cos(angle),np.sin(angle)])
        states.append(analytic_adiabatic_vectors(R)[:,0])

    graph=ElectronicGaugeGraph(1)
    for i in range(n):
        j=(i+1)%n
        overlap=np.vdot(states[i],states[j])
        graph.add_overlap(i,j,np.array([[overlap]]))

    W=graph.wilson_loop(list(range(n)))
    assert np.real(W[0,0]) < -0.999999
    assert abs(np.imag(W[0,0])) < 1e-12
