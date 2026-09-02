from dataclasses import dataclass
import numpy as np
from .gauge_graph import ElectronicGaugeGraph
from .graph_electronic import GraphElectronicRegistry
from .ci2d import LVC2DParameters,adiabatic_energies_2d,adiabatic_gradients_2d,analytic_adiabatic_vectors,vector_nac_2d

@dataclass
class ElectronicFramePoint:
    q: np.ndarray; frame: np.ndarray; energies: np.ndarray; gradients: np.ndarray; nac: np.ndarray; mass_matrix: np.ndarray
    def __post_init__(self):
        self.q=np.asarray(self.q,float); self.frame=np.asarray(self.frame,complex); self.energies=np.asarray(self.energies,float)
        self.gradients=np.asarray(self.gradients,float); self.nac=np.asarray(self.nac,float); self.mass_matrix=np.asarray(self.mass_matrix,float)
        ns=len(self.energies); nq=len(self.q)
        if self.frame.shape!=(ns,ns): raise ValueError('frame must have shape (nstate,nstate)')
        if self.gradients.shape!=(ns,nq): raise ValueError('gradients must have shape (nstate,nq)')
        if self.nac.shape!=(ns,ns,nq): raise ValueError('nac must have shape (nstate,nstate,nq)')
        if self.mass_matrix.shape!=(nq,nq): raise ValueError('mass_matrix must have shape (nq,nq)')
        if not np.allclose(self.frame.conj().T@self.frame,np.eye(ns),atol=1e-10): raise ValueError('frame must be unitary')
        if not np.allclose(self.mass_matrix,self.mass_matrix.T,atol=1e-12): raise ValueError('mass_matrix must be symmetric')

class AnalyticCI2DFrameProvider:
    def __init__(self,nuclear_mass_au=20.0,params=LVC2DParameters()): self.nuclear_mass_au=float(nuclear_mass_au); self.params=params
    def evaluate(self,q):
        q=np.asarray(q,float)
        if q.shape!=(2,): raise ValueError('AnalyticCI2DFrameProvider expects q=(x,y)')
        return ElectronicFramePoint(q.copy(),analytic_adiabatic_vectors(q,self.params),adiabatic_energies_2d(q[0],q[1],self.params),adiabatic_gradients_2d(q,self.params),vector_nac_2d(q,self.params),self.nuclear_mass_au*np.eye(2))

class IncrementalElectronicGraph:
    def __init__(self,dimension):
        self.graph=ElectronicGaugeGraph(dimension); self.registry=GraphElectronicRegistry(self.graph); self.frames={}; self.coordinates={}; self.mass_matrices={}
    @property
    def dimension(self): return self.graph.dimension
    def add_frame_point(self,node,point,connect_to=()):
        if node in self.frames: raise ValueError(f'node {node!r} already exists')
        if point.frame.shape!=(self.dimension,self.dimension): raise ValueError('electronic dimension differs from graph dimension')
        self.graph.add_node(node); self.registry.add_adiabatic_data(node,point.energies,point.gradients,point.nac)
        self.frames[node]=point.frame.copy(); self.coordinates[node]=point.q.copy(); self.mass_matrices[node]=point.mass_matrix.copy()
        for other in connect_to: self.connect(node,other)
    def connect(self,u,v,weight=None):
        if u not in self.frames or v not in self.frames: raise KeyError('both nodes must exist before connecting')
        overlap=self.frames[u].conj().T@self.frames[v]; self.graph.add_overlap(u,v,overlap,weight=weight); return overlap
    def add_from_provider(self,node,q,provider,connect_to=()):
        point=provider.evaluate(q); self.add_frame_point(node,point,connect_to); return point
    def temporal_link(self,old_node,new_node): return self.graph.link(old_node,new_node)
    def summary(self): return {'nodes':len(self.graph.nodes),'edges':len(self.graph.edges()),'fundamental_cycles':max(0,len(self.graph.edges())-len(self.graph.nodes)+1)}
