from dataclasses import dataclass
import numpy as np


def adiabatic_hamiltonian_matrix(energies):
    E = np.asarray(energies, dtype=float)
    if E.ndim != 1:
        raise ValueError("energies must be one-dimensional")
    return np.diag(E).astype(complex)


def derivative_hamiltonian_matrices(energies, gradients, nac):
    """Construct F_alpha = <phi_i|dH/dq_alpha|phi_j>.

    In an adiabatic basis,

        F_ii,alpha = dE_i/dq_alpha
        F_ij,alpha = (E_j-E_i) d_ij,alpha  for i != j.

    Output shape: (nq,nstate,nstate).
    """
    E = np.asarray(energies, dtype=float)
    G = np.asarray(gradients, dtype=float)
    D = np.asarray(nac, dtype=float)

    ns = len(E)
    if G.ndim != 2 or G.shape[0] != ns:
        raise ValueError("gradients must have shape (nstate,nq)")
    nq = G.shape[1]
    if D.shape != (ns, ns, nq):
        raise ValueError("nac must have shape (nstate,nstate,nq)")

    F = np.zeros((nq, ns, ns), dtype=complex)
    for a in range(ns):
        F[:, a, a] = G[a]
    for i in range(ns):
        for j in range(ns):
            if i == j:
                continue
            F[:, i, j] = (E[j] - E[i]) * D[i, j]
    return F


def rotate_operator(operator, gauge):
    A = np.asarray(operator, dtype=complex)
    G = np.asarray(gauge, dtype=complex)
    return G.conj().T @ A @ G


def rotate_operator_field(field, gauge):
    F = np.asarray(field, dtype=complex)
    return np.asarray([rotate_operator(A, gauge) for A in F])


def rotate_coefficients(coefficients, gauge):
    """If Phi' = Phi G, then c' = G^dagger c."""
    c = np.asarray(coefficients, dtype=complex)
    G = np.asarray(gauge, dtype=complex)
    return G.conj().T @ c


@dataclass
class ElectronicOperatorNode:
    hamiltonian: np.ndarray
    derivative_hamiltonians: np.ndarray

    def __post_init__(self):
        self.hamiltonian = np.asarray(self.hamiltonian, dtype=complex)
        self.derivative_hamiltonians = np.asarray(
            self.derivative_hamiltonians, dtype=complex
        )
        ns = self.hamiltonian.shape[0]
        if self.hamiltonian.shape != (ns, ns):
            raise ValueError("hamiltonian must be square")
        if self.derivative_hamiltonians.ndim != 3:
            raise ValueError("derivative_hamiltonians must have shape (nq,ns,ns)")
        if self.derivative_hamiltonians.shape[1:] != (ns, ns):
            raise ValueError("derivative_hamiltonians have incompatible shape")


class GraphElectronicRegistry:
    """Electronic operators attached to the nodes of an ElectronicGaugeGraph."""

    def __init__(self, graph):
        self.graph = graph
        self.data = {}

    def add_adiabatic_data(self, node, energies, gradients, nac):
        H = adiabatic_hamiltonian_matrix(energies)
        F = derivative_hamiltonian_matrices(energies, gradients, nac)
        self.add_operator_data(node, H, F)

    def add_operator_data(self, node, hamiltonian, derivative_hamiltonians):
        if node not in self.graph.nodes:
            self.graph.add_node(node)
        data = ElectronicOperatorNode(hamiltonian, derivative_hamiltonians)
        if data.hamiltonian.shape[0] != self.graph.dimension:
            raise ValueError("operator dimension differs from graph dimension")
        self.data[node] = data

    def transport_coefficients(self, source, target, coefficients, path=None):
        return self.graph.transport_coefficients(source, target, coefficients, path=path)

    def pair_factors(self, node_i, coeff_i, node_j, coeff_j, reference_node):
        """Gauge-covariant scalar electronic factors in a chosen reference frame."""
        if reference_node not in self.data:
            raise KeyError("reference node has no operator data")
        vi = self.transport_coefficients(node_i, reference_node, coeff_i)
        vj = self.transport_coefficients(node_j, reference_node, coeff_j)
        op = self.data[reference_node]

        overlap = np.vdot(vi, vj)
        potential = np.vdot(vi, op.hamiltonian @ vj)
        derivatives = np.asarray([
            np.vdot(vi, F @ vj)
            for F in op.derivative_hamiltonians
        ])
        return {
            "overlap": overlap,
            "potential": potential,
            "derivative_hamiltonian": derivatives,
            "state_i_at_reference": vi,
            "state_j_at_reference": vj,
        }

    def gauge_transformed(self, gauges):
        """Return a covariantly transformed registry and graph links.

        The transformed graph is rebuilt from transformed unitary links.  Operator
        matrices transform as G^dagger A G.
        """
        from .gauge_graph import ElectronicGaugeGraph

        new_graph = ElectronicGaugeGraph(self.graph.dimension)
        for node in self.graph.nodes:
            new_graph.add_node(node)

        for edge in self.graph.edges():
            Gu = np.asarray(gauges[edge.u], dtype=complex)
            Gv = np.asarray(gauges[edge.v], dtype=complex)
            O = Gu.conj().T @ self.graph.overlap(edge.u, edge.v) @ Gv
            new_graph.add_overlap(edge.u, edge.v, O, weight=edge.weight)

        out = GraphElectronicRegistry(new_graph)
        for node, data in self.data.items():
            G = np.asarray(gauges[node], dtype=complex)
            out.add_operator_data(
                node,
                rotate_operator(data.hamiltonian, G),
                rotate_operator_field(data.derivative_hamiltonians, G),
            )
        return out
