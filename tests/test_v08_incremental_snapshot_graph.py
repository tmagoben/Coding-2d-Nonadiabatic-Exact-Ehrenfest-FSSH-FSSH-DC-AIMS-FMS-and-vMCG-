import numpy as np
from gaussian_dynamics.incremental_snapshot_graph import IncrementalSnapshotGaugeGraph
class Snap:
    def __init__(self,label,nroots=2):self.label=label;self.nroots=nroots

def test_incremental_snapshot_graph_uses_overlap_engine():
    frames={'a':np.eye(2,dtype=complex),'b':np.array([[.8,-.6],[.6,.8]],complex),'c':np.array([[.6,-.8],[.8,.6]],complex)};sn={k:Snap(k) for k in frames};engine=lambda old,new:frames[old.label].conj().T@frames[new.label];b=IncrementalSnapshotGaugeGraph(2,engine);E=np.array([.1,.3]);G=np.array([[.01,0],[-.02,0]]);D=np.zeros((2,2,2));D[0,1,0]=.2;D[1,0,0]=-.2
    b.add_adiabatic_node('a',sn['a'],E,G,D);b.add_adiabatic_node('b',sn['b'],E,G,D,['a']);b.add_adiabatic_node('c',sn['c'],E,G,D,['a','b']);assert len(b.graph.edges())==3;assert np.allclose(b.graph.wilson_loop(['a','b','c']),np.eye(2),atol=1e-12)

def test_edge_diagnostics_report_nonunitary_raw_overlap():
    O=np.array([[.95,0],[0,.8]],complex);b=IncrementalSnapshotGaugeGraph(2,lambda a,c:O);H=np.diag([.1,.2]);F=np.zeros((1,2,2));b.add_node(0,Snap('0'),H,F);b.add_node(1,Snap('1'),H,F,[0]);r=b.edge_diagnostics()[0];assert np.allclose(r['singular_values'],[.95,.8]) and r['unitarity_defect']>0
