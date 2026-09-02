import numpy as np

from gaussian_dynamics.dynamic_graph_aims import DynamicGraphTBF
from gaussian_dynamics.locality_graph_v16 import (
    conservative_position_overlap_bound,
    LocalityGraphSettings,
    PersistentGaussianLocalityGraph,
)
from gaussian_dynamics.pair_cache_v15 import GaussianPairCache


def _tbf(uid,q,p,A):
    return DynamicGraphTBF(
        uid=uid,state=0,
        q=np.asarray(q,float),
        p=np.asarray(p,float),
        A=np.asarray(A,float),
        node=("t",uid),
    )


def test_conservative_position_bound_exceeds_exact_overlap():
    rng=np.random.default_rng(1234)

    for k in range(30):
        Q,_=np.linalg.qr(rng.normal(size=(2,2)))
        vals=np.diag(rng.uniform(0.5,2.0,size=2))
        Ai=Q@vals@Q.T

        Q,_=np.linalg.qr(rng.normal(size=(2,2)))
        vals=np.diag(rng.uniform(0.5,2.0,size=2))
        Aj=Q@vals@Q.T

        a=_tbf(
            2*k,
            rng.normal(size=2),
            rng.normal(size=2),
            Ai,
        )
        b=_tbf(
            2*k+1,
            rng.normal(size=2),
            rng.normal(size=2),
            Aj,
        )

        bound=conservative_position_overlap_bound(a,b)
        exact=abs(GaussianPairCache([a,b]).pair(0,1).overlap)

        assert exact<=bound+2e-13


def test_locality_graph_screens_far_pair_without_exact_pair_solve():
    basis=[
        _tbf(0,[0,0],[0,0],np.eye(2)),
        _tbf(1,[8,0],[0,0],np.eye(2)),
    ]
    graph=PersistentGaussianLocalityGraph(
        LocalityGraphSettings(
            enter_overlap=1e-4,
            exit_overlap=5e-5,
        )
    )
    update=graph.update(basis)

    assert update.active_edges==()
    assert update.screened_pairs==1
    assert update.exact_pair_checks==0
    assert update.cache.stats.canonical_solves==0


def test_edge_hysteresis_uses_smaller_exit_threshold():
    # Equal unit widths: |S| = exp(-dq^2/4) at zero momentum.
    enter=0.20
    exit_=0.10
    graph=PersistentGaussianLocalityGraph(
        LocalityGraphSettings(
            enter_overlap=enter,
            exit_overlap=exit_,
        )
    )

    # First geometry: overlap > enter -> edge appears.
    b0=[
        _tbf(0,[0,0],[0,0],np.eye(2)),
        _tbf(1,[2.0,0],[0,0],np.eye(2)),
    ]
    u0=graph.update(b0)
    assert len(u0.active_edges)==1

    # Move apart enough that overlap falls below enter but remains above exit.
    b1=[
        _tbf(0,[0,0],[0,0],np.eye(2)),
        _tbf(1,[2.7,0],[0,0],np.eye(2)),
    ]
    u1=graph.update(b1)
    assert len(u1.active_edges)==1
    assert u1.retained_edges==1

    # Move farther: overlap below exit -> edge removed.
    b2=[
        _tbf(0,[0,0],[0,0],np.eye(2)),
        _tbf(1,[3.2,0],[0,0],np.eye(2)),
    ]
    u2=graph.update(b2)
    assert len(u2.active_edges)==0
    assert u2.exited_edges==1
