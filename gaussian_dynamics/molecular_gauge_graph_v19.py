from dataclasses import dataclass
import numpy as np

from .gauge_graph import ElectronicGaugeGraph
from .graph_electronic import GraphElectronicRegistry
from .graph_gaussian import (
    GraphGaussianTBF,
    build_static_graph_gaussian_matrices,
)


@dataclass
class MolecularCentroidGraphV19:
    registry: GraphElectronicRegistry
    basis: list
    pair_centroid_nodes: dict
    center_nodes: tuple
    mass_matrix: np.ndarray
    provider_nodes: dict

    def reference_selector(self,i,j):
        if i==j:
            return self.center_nodes[i]
        key=(min(int(i),int(j)),max(int(i),int(j)))
        return self.pair_centroid_nodes[key]

    def matrices(self):
        return build_static_graph_gaussian_matrices(
            self.basis,
            self.registry,
            self.mass_matrix,
            self.reference_selector,
        )


def build_molecular_centroid_graph_v19(
    local_basis,
    provider,
):
    """Build a center-centroid electronic gauge graph for a molecular Gaussian basis.

    Each local TBF is converted to a GraphGaussianTBF carrying the unit vector
    corresponding to its tracked local adiabatic state. Every pair centroid becomes an
    electronic-operator node connected to both TBF centers using cross-geometry
    electronic overlaps.

    For a real PySCF backend the provider should expose overlap-capable tracked
    snapshots. The deterministic v0.19 backend uses finite state vectors.
    """
    local_basis=list(local_basis)
    if not local_basis:
        raise ValueError("local_basis cannot be empty.")

    center_nodes=tuple(
        f"center-{i}"
        for i in range(len(local_basis))
    )
    center_snapshots={}
    ns=None
    mass=None

    for i,b in enumerate(local_basis):
        snap=provider.evaluate_snapshot(b.q)
        center_snapshots[center_nodes[i]]=snap
        if ns is None:
            ns=len(snap.point.energies)
            mass=snap.point.mass_matrix_q_au.copy()
        elif len(snap.point.energies)!=ns:
            raise ValueError("number of electronic states changed across centers.")
        if not np.allclose(
            snap.point.mass_matrix_q_au,
            mass,atol=1e-10,rtol=1e-10
        ):
            raise ValueError(
                "v0.19 centroid graph currently requires one constant generalized mass matrix."
            )

    graph=ElectronicGaugeGraph(ns)
    registry=GraphElectronicRegistry(graph)
    provider_nodes=dict(center_snapshots)

    for node,snap in center_snapshots.items():
        registry.add_adiabatic_data(
            node,
            snap.point.energies,
            snap.point.gradients_q,
            snap.point.nac_q,
        )

    pair_centroids={}
    for i in range(len(local_basis)):
        for j in range(i+1,len(local_basis)):
            node=f"centroid-{i}-{j}"
            qbar=0.5*(
                np.asarray(local_basis[i].q,float)
                +np.asarray(local_basis[j].q,float)
            )
            snap=provider.evaluate_snapshot(qbar)
            pair_centroids[(i,j)]=node
            provider_nodes[node]=snap
            registry.add_adiabatic_data(
                node,
                snap.point.energies,
                snap.point.gradients_q,
                snap.point.nac_q,
            )

    def overlap(a,b):
        A=provider_nodes[a].state_vectors
        B=provider_nodes[b].state_vectors
        if A is None or B is None:
            raise ValueError(
                "Finite state_vectors are required by this v0.19 validation graph. "
                "For PySCF use the snapshot overlap engine bridge."
            )
        return A.conj().T@B

    for (i,j),centroid in pair_centroids.items():
        graph.add_overlap(
            center_nodes[i],centroid,
            overlap(center_nodes[i],centroid),
        )
        graph.add_overlap(
            center_nodes[j],centroid,
            overlap(center_nodes[j],centroid),
        )

    # Ensure a one-TBF graph still has the center node.
    for node in center_nodes:
        graph.add_node(node)

    gbasis=[]
    for i,b in enumerate(local_basis):
        c=np.zeros(ns,dtype=complex)
        c[int(b.state)]=1.0
        gbasis.append(
            GraphGaussianTBF(
                node=center_nodes[i],
                q=np.asarray(b.q,float).copy(),
                p=np.asarray(b.p,float).copy(),
                A=np.asarray(b.A,float).copy(),
                electronic_coefficients=c,
            )
        )

    return MolecularCentroidGraphV19(
        registry=registry,
        basis=gbasis,
        pair_centroid_nodes=pair_centroids,
        center_nodes=center_nodes,
        mass_matrix=mass,
        provider_nodes=provider_nodes,
    )
