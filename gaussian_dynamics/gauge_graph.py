from dataclasses import dataclass
from collections import deque
import numpy as np

from .finite_manifold_transport_v233 import (
    CONSUMER_OVERLAP_POLICY_V233,
    certified_transport_from_overlap_v233,
)


@dataclass(frozen=True)
class GaugeEdge:
    """One undirected electronic-overlap edge.

    The stored orientation is u -> v with

        O_uv = <Phi_u | Phi_v>

    and link_uv is the nearest unitary (polar factor) to O_uv.
    The reverse link is link_uv^dagger.
    """
    u: object
    v: object
    overlap_uv: np.ndarray
    link_uv: np.ndarray
    weight: float


def nearest_unitary(matrix):
    """Generic unitary polar projection used by graph optimization.

    Spectral synchronization and coordinate refinement apply this function to
    sums/eigenvector blocks that are not physical cross-geometry overlaps and
    therefore are not contractions. Physical edge overlaps are certified
    separately in :meth:`ElectronicGaugeGraph.add_overlap`.
    """
    matrix = np.asarray(matrix, dtype=complex)
    if (
        matrix.ndim != 2
        or matrix.shape[0] < 1
        or matrix.shape[0] != matrix.shape[1]
        or not np.all(np.isfinite(matrix))
    ):
        raise ValueError("polar projection requires a finite nonempty square matrix")
    left, _, right_h = np.linalg.svd(matrix, full_matrices=False)
    return left @ right_h


class ElectronicGaugeGraph:
    """Discrete electronic gauge connection on a graph.

    Local electronic frames live at graph nodes.  Edge links transform as

        U_uv -> G_u^dagger U_uv G_v.

    Wilson-loop spectra are therefore gauge invariant.
    """

    def __init__(self, dimension):
        self.dimension = int(dimension)
        if self.dimension <= 0:
            raise ValueError("dimension must be positive")
        self.nodes = set()
        self._edges = {}
        self._adj = {}

    def add_node(self, node):
        self.nodes.add(node)
        self._adj.setdefault(node, set())

    def add_overlap(self, u, v, overlap, weight=None):
        if u == v:
            raise ValueError("self edges are not needed")
        O = np.asarray(overlap, dtype=complex)
        if O.shape != (self.dimension, self.dimension):
            raise ValueError("overlap has incompatible dimension")

        self.add_node(u)
        self.add_node(v)

        link = certified_transport_from_overlap_v233(
            O, policy=CONSUMER_OVERLAP_POLICY_V233
        ).right_to_left_transport.copy()
        if weight is None:
            # Mean singular value is a simple edge-confidence diagnostic.
            weight = float(np.mean(np.linalg.svd(O, compute_uv=False)))
        weight = float(weight)
        if weight <= 0.0:
            raise ValueError("edge weight must be positive")

        key = frozenset((u, v))
        self._edges[key] = GaugeEdge(u, v, O.copy(), link, weight)
        self._adj[u].add(v)
        self._adj[v].add(u)

    def neighbors(self, node):
        return tuple(self._adj[node])

    def edges(self):
        return tuple(self._edges.values())

    def _edge(self, u, v):
        key = frozenset((u, v))
        if key not in self._edges:
            raise KeyError(f"No edge between {u!r} and {v!r}")
        return self._edges[key]

    def overlap(self, u, v):
        edge = self._edge(u, v)
        if edge.u == u and edge.v == v:
            return edge.overlap_uv.copy()
        return edge.overlap_uv.conj().T.copy()

    def link(self, u, v):
        """Return U_uv = nearest-unitary(<Phi_u|Phi_v>)."""
        edge = self._edge(u, v)
        if edge.u == u and edge.v == v:
            return edge.link_uv.copy()
        return edge.link_uv.conj().T.copy()

    def weight(self, u, v):
        return self._edge(u, v).weight

    def path(self, source, target):
        """Shortest unweighted graph path by BFS."""
        if source not in self.nodes or target not in self.nodes:
            raise KeyError("source/target node missing")
        if source == target:
            return [source]

        parent = {source: None}
        q = deque([source])
        while q:
            u = q.popleft()
            for v in self._adj[u]:
                if v in parent:
                    continue
                parent[v] = u
                if v == target:
                    q.clear()
                    break
                q.append(v)

        if target not in parent:
            raise ValueError("graph is disconnected")

        out = [target]
        while out[-1] != source:
            out.append(parent[out[-1]])
        out.reverse()
        return out

    def transport_matrix(self, source, target, path=None):
        """Map electronic coefficients from source-frame coordinates to target.

        If U_uv = <Phi_u|Phi_v>, then coefficients transform from v to u as

            c_u = U_uv c_v.

        Along source=n0 -> ... -> nk=target, the required product is

            U_{nk,n{k-1}} ... U_{n1,n0}.
        """
        if path is None:
            path = self.path(source, target)
        if path[0] != source or path[-1] != target:
            raise ValueError("path endpoints do not match source/target")

        T = np.eye(self.dimension, dtype=complex)
        for a, b in zip(path[:-1], path[1:]):
            T = self.link(b, a) @ T
        return T

    def transport_coefficients(self, source, target, coefficients, path=None):
        c = np.asarray(coefficients, dtype=complex)
        if c.shape != (self.dimension,):
            raise ValueError("coefficient vector has incompatible dimension")
        return self.transport_matrix(source, target, path=path) @ c

    def wilson_loop(self, cycle):
        """Ordered Wilson product around a closed node cycle.

        `cycle` may omit the repeated first node.  For n0,n1,...,nk,

            W = U_n0,n1 U_n1,n2 ... U_nk,n0.

        Under local gauges W -> G_n0^dagger W G_n0, so its eigenvalues and trace
        are gauge invariant.
        """
        cycle = list(cycle)
        if len(cycle) < 2:
            raise ValueError("cycle must contain at least two distinct nodes")
        if cycle[0] == cycle[-1]:
            cycle = cycle[:-1]

        W = np.eye(self.dimension, dtype=complex)
        closed = cycle + [cycle[0]]
        for u, v in zip(closed[:-1], closed[1:]):
            W = W @ self.link(u, v)
        return W

    def spanning_tree(self, root):
        if root not in self.nodes:
            raise KeyError(root)
        parent = {root: None}
        q = deque([root])
        while q:
            u = q.popleft()
            for v in sorted(self._adj[u], key=repr):
                if v in parent:
                    continue
                parent[v] = u
                q.append(v)
        if len(parent) != len(self.nodes):
            raise ValueError("graph is disconnected")
        return parent

    def spanning_tree_gauges(self, root):
        """Choose node gauges that make every spanning-tree link the identity."""
        parent = self.spanning_tree(root)
        gauges = {root: np.eye(self.dimension, dtype=complex)}

        pending = [n for n in parent if n != root]
        while pending:
            progress = False
            for node in list(pending):
                p = parent[node]
                if p not in gauges:
                    continue
                # G_p^dag U_p,node G_node = I
                gauges[node] = self.link(p, node).conj().T @ gauges[p]
                pending.remove(node)
                progress = True
            if not progress:
                raise RuntimeError("failed to construct tree gauges")
        return gauges, parent

    def transformed_link(self, u, v, gauges):
        Gu = np.asarray(gauges[u], dtype=complex)
        Gv = np.asarray(gauges[v], dtype=complex)
        return Gu.conj().T @ self.link(u, v) @ Gv

    def gauge_objective(self, gauges):
        I = np.eye(self.dimension, dtype=complex)
        total = 0.0
        for edge in self.edges():
            L = self.transformed_link(edge.u, edge.v, gauges)
            total += edge.weight * np.linalg.norm(L - I, ord="fro") ** 2
        return float(total)

    def spectral_synchronize(self, root):
        """Spectral relaxation for unitary synchronization.

        Build the block Hermitian connection matrix

            C_uv = w_uv U_uv,
            C_vu = w_uv U_uv^dagger.

        The top `dimension` eigenvectors give a relaxed solution.  Each node block is
        projected back to U(m) by its polar factor, then the global right-unitary
        freedom is fixed by anchoring the requested root to I.
        """
        if root not in self.nodes:
            raise KeyError(root)

        nodes = sorted(self.nodes, key=repr)
        index = {node: i for i, node in enumerate(nodes)}
        m = self.dimension
        C = np.zeros((len(nodes)*m, len(nodes)*m), dtype=complex)

        for edge in self.edges():
            i = index[edge.u]
            j = index[edge.v]
            si = slice(i*m, (i+1)*m)
            sj = slice(j*m, (j+1)*m)
            U = self.link(edge.u, edge.v)
            C[si, sj] += edge.weight * U
            C[sj, si] += edge.weight * U.conj().T

        evals, evecs = np.linalg.eigh(C)
        Y = evecs[:, -m:]

        gauges = {}
        for node in nodes:
            i = index[node]
            block = Y[i*m:(i+1)*m, :]
            gauges[node] = nearest_unitary(block)

        anchor = gauges[root].conj().T
        gauges = {node: G @ anchor for node, G in gauges.items()}
        gauges[root] = np.eye(m, dtype=complex)
        return gauges

    def _coordinate_refine(self, gauges, root, max_iter, tolerance):
        gauges = {node: np.asarray(G, dtype=complex).copy() for node, G in gauges.items()}
        I = np.eye(self.dimension, dtype=complex)
        gauges[root] = I

        for _ in range(int(max_iter)):
            max_change = 0.0
            for u in sorted(self.nodes, key=repr):
                if u == root:
                    continue
                A = np.zeros((self.dimension, self.dimension), dtype=complex)
                for v in self._adj[u]:
                    A += self.weight(u, v) * self.link(u, v) @ gauges[v]
                if np.linalg.norm(A) == 0.0:
                    continue
                new = nearest_unitary(A)
                max_change = max(
                    max_change,
                    np.linalg.norm(new - gauges[u], ord="fro"),
                )
                gauges[u] = new
            gauges[root] = I
            if max_change < tolerance:
                break
        return gauges

    def synchronize(
        self,
        root,
        max_iter=100,
        tolerance=1e-10,
        restarts=2,
        seed=0,
    ):
        """Anchored unitary synchronization with deterministic multi-start refinement.

        The objective is

            sum_(u,v) w_uv ||G_u^dag U_uv G_v - I||_F^2.

        A spectral and a spanning-tree initialization are supplemented by a small
        number of seeded random-unitary starts.  Each is refined by transparent block
        coordinate descent and the lowest-objective result is returned.

        Exact zero is possible only for a flat/holonomy-free graph connection.
        Nontrivial Wilson loops leave an irreducible residual.
        """
        if root not in self.nodes:
            raise KeyError(root)

        candidates = []
        spectral = self.spectral_synchronize(root)
        candidates.append(spectral)

        tree, _ = self.spanning_tree_gauges(root)
        candidates.append(tree)

        rng = np.random.default_rng(seed)
        m = self.dimension
        for _ in range(int(restarts)):
            trial = {root: np.eye(m, dtype=complex)}
            for node in self.nodes:
                if node == root:
                    continue
                X = rng.normal(size=(m, m)) + 1j*rng.normal(size=(m, m))
                Q, _ = np.linalg.qr(X)
                trial[node] = Q
            candidates.append(trial)

        best = None
        best_obj = np.inf
        for trial in candidates:
            refined = self._coordinate_refine(
                trial,
                root=root,
                max_iter=max_iter,
                tolerance=tolerance,
            )
            obj = self.gauge_objective(refined)
            if obj < best_obj:
                best_obj = obj
                best = refined
        return best

    def fundamental_cycles(self, root):
        """Fundamental cycle basis generated by chords relative to a BFS tree."""
        parent = self.spanning_tree(root)
        tree_edges = {
            frozenset((node, p))
            for node, p in parent.items()
            if p is not None
        }

        tree_adj = {n: set() for n in self.nodes}
        for edge in tree_edges:
            u, v = tuple(edge)
            tree_adj[u].add(v)
            tree_adj[v].add(u)

        def tree_path(a, b):
            par = {a: None}
            q = deque([a])
            while q:
                x = q.popleft()
                if x == b:
                    break
                for y in tree_adj[x]:
                    if y not in par:
                        par[y] = x
                        q.append(y)
            out = [b]
            while out[-1] != a:
                out.append(par[out[-1]])
            out.reverse()
            return out

        cycles = []
        for edge in self.edges():
            key = frozenset((edge.u, edge.v))
            if key in tree_edges:
                continue
            path = tree_path(edge.u, edge.v)
            cycles.append(path)
        return tuple(cycles)
