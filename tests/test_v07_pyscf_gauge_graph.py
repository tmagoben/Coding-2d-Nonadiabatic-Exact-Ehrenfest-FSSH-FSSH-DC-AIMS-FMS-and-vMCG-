import numpy as np

from gaussian_dynamics.pyscf_gauge_graph import (
    build_snapshot_gauge_graph,
    edge_overlap_diagnostics,
    tbf_centroid_edge_pairs,
)


class FakeSnapshot:
    def __init__(self,nroots,label):
        self.nroots=nroots
        self.label=label


def test_snapshot_graph_builder_and_diagnostics():
    snaps={k:FakeSnapshot(2,k) for k in ("a","b","c")}
    frames={
        "a":np.eye(2),
        "b":np.array([[0.8,-0.6],[0.6,0.8]]),
        "c":np.array([[0.6,-0.8],[0.8,0.6]]),
    }

    def engine(s1,s2):
        return frames[s1.label].T @ frames[s2.label]

    graph=build_snapshot_gauge_graph(
        snaps,
        [("a","b"),("b","c"),("c","a")],
        overlap_engine=engine,
    )

    assert graph.dimension==2
    assert len(graph.edges())==3
    assert np.allclose(graph.wilson_loop(["a","b","c"]),np.eye(2),atol=1e-12)

    diag=edge_overlap_diagnostics(graph)
    assert len(diag)==3
    assert max(d["unitarity_defect"] for d in diag) < 1e-12


def test_tbf_centroid_connectivity_creates_expected_edges():
    edges=tbf_centroid_edge_pairs(
        ["t0","t1","t2"],
        {(0,1):"c01",(1,2):"c12",(0,2):"c02"},
    )
    assert len(edges)==6
    assert ("t0","c01") in edges
    assert ("t1","c01") in edges
    assert ("t2","c02") in edges
