import numpy as np

from .gauge_graph import ElectronicGaugeGraph
from .pyscf_wavefunction_overlap import casscf_state_overlap_matrix


def build_snapshot_gauge_graph(
    snapshots,
    edge_pairs,
    overlap_engine=None,
):
    """Build a gauge graph from PySCF/CASSCF wavefunction snapshots.

    Parameters
    ----------
    snapshots
        Mapping node_id -> CASSCFWavefunctionSnapshot.
    edge_pairs
        Iterable of (u,v) node pairs whose many-electron overlap should be evaluated.
    overlap_engine
        Defaults to `casscf_state_overlap_matrix`.  Injectable for tests.
    """
    snapshots = dict(snapshots)
    if not snapshots:
        raise ValueError("snapshots must not be empty")

    first = next(iter(snapshots.values()))
    dimension = int(first.nroots)
    for node, snap in snapshots.items():
        if int(snap.nroots) != dimension:
            raise ValueError(f"snapshot {node!r} has a different number of roots")

    engine = overlap_engine or casscf_state_overlap_matrix
    graph = ElectronicGaugeGraph(dimension)
    for node in snapshots:
        graph.add_node(node)

    for u, v in edge_pairs:
        O = np.asarray(engine(snapshots[u], snapshots[v]), dtype=complex)
        graph.add_overlap(u, v, O)

    return graph


def edge_overlap_diagnostics(graph):
    """Return singular-value and unitarity-defect diagnostics for every edge."""
    out = []
    I = np.eye(graph.dimension, dtype=complex)
    for edge in graph.edges():
        O = graph.overlap(edge.u, edge.v)
        singular_values = np.linalg.svd(O, compute_uv=False)
        defect = np.linalg.norm(O.conj().T @ O - I, ord="fro")
        out.append({
            "u": edge.u,
            "v": edge.v,
            "weight": edge.weight,
            "singular_values": singular_values,
            "unitarity_defect": float(defect),
        })
    return out


def tbf_centroid_edge_pairs(tbf_nodes, pair_centroid_nodes):
    """Connectivity helper for branched Gaussian bases.

    pair_centroid_nodes maps unordered TBF-index pairs to centroid node IDs, e.g.

        {(0,1): "c01", (0,2): "c02", (1,2): "c12"}.

    The returned graph connects each centroid to both TBF-center nodes.  These
    center-centroid-center paths generate loops whenever multiple TBF pairs are present.
    """
    tbf_nodes = list(tbf_nodes)
    edges = []
    for pair, centroid in pair_centroid_nodes.items():
        if len(pair) != 2:
            raise ValueError("pair keys must contain two TBF indices")
        i, j = pair
        edges.append((tbf_nodes[i], centroid))
        edges.append((tbf_nodes[j], centroid))
    return tuple(edges)
