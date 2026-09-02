import numpy as np
from .gauge_graph import ElectronicGaugeGraph
from .graph_electronic import GraphElectronicRegistry,adiabatic_hamiltonian_matrix,derivative_hamiltonian_matrices
from .pyscf_wavefunction_overlap import casscf_state_overlap_matrix

class IncrementalSnapshotGaugeGraph:
    def __init__(self,dimension,overlap_engine=None):
        self.graph=ElectronicGaugeGraph(int(dimension)); self.registry=GraphElectronicRegistry(self.graph); self.snapshots={}; self.overlap_engine=overlap_engine or casscf_state_overlap_matrix
    def add_node(self,node,snapshot,hamiltonian,derivative_hamiltonians,connect_to=()):
        if node in self.snapshots: raise ValueError(f'node {node!r} already exists')
        if int(snapshot.nroots)!=self.graph.dimension: raise ValueError('snapshot root count differs from graph dimension')
        self.graph.add_node(node); self.registry.add_operator_data(node,hamiltonian,derivative_hamiltonians); self.snapshots[node]=snapshot; out={}
        for other in connect_to:
            if other not in self.snapshots: raise KeyError(other)
            O=np.asarray(self.overlap_engine(self.snapshots[other],snapshot),complex); self.graph.add_overlap(other,node,O);out[other]=O
        return out
    def add_adiabatic_node(self,node,snapshot,energies,gradients,nac,connect_to=()):
        return self.add_node(node,snapshot,adiabatic_hamiltonian_matrix(energies),derivative_hamiltonian_matrices(energies,gradients,nac),connect_to)
    def add_cartesian_point(self,node,snapshot,cartesian_point,connect_to=()):
        p=cartesian_point.validate(); ns=len(p.energies); return self.add_adiabatic_node(node,snapshot,p.energies,p.gradients_cart.reshape(ns,-1),p.nac_cart.reshape(ns,ns,-1),connect_to)
    def edge_diagnostics(self):
        I=np.eye(self.graph.dimension,dtype=complex); rows=[]
        for e in self.graph.edges():
            O=self.graph.overlap(e.u,e.v); rows.append({'u':e.u,'v':e.v,'singular_values':np.linalg.svd(O,compute_uv=False),'unitarity_defect':float(np.linalg.norm(O.conj().T@O-I,ord='fro')),'weight':e.weight})
        return rows
