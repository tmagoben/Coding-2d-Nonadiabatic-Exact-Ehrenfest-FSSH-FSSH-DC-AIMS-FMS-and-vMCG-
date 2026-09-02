import numpy as np
from gaussian_dynamics.dynamic_gauge_graph import AnalyticCI2DFrameProvider,IncrementalElectronicGraph

def test_incremental_graph_edges_form_cycle():
    p=AnalyticCI2DFrameProvider();m=IncrementalElectronicGraph(2);m.add_from_provider('a',np.array([-.8,.6]),p);m.add_from_provider('b',np.array([.7,.5]),p);m.add_from_provider('c',np.array([-.05,.55]),p,connect_to=['a','b']);m.connect('a','b')
    s=m.summary();assert s['nodes']==3 and s['edges']==3 and s['fundamental_cycles']==1
    W=m.graph.wilson_loop(['a','c','b']);assert np.allclose(W.conj().T@W,np.eye(2),atol=1e-12)

def test_temporal_link_is_unitary():
    p=AnalyticCI2DFrameProvider();m=IncrementalElectronicGraph(2);m.add_from_provider('old',np.array([-.8,.6]),p);m.add_from_provider('new',np.array([-.79,.605]),p,connect_to=['old']);U=m.temporal_link('old','new');assert np.allclose(U.conj().T@U,np.eye(2),atol=1e-12)
