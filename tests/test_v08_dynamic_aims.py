import numpy as np
from gaussian_dynamics.dynamic_graph_aims import DynamicGraphTBF,run_dynamic_graph_aims

def test_dynamic_graph_aims_grows_spawns_and_preserves_norm():
    b=DynamicGraphTBF(0,1,np.array([.55,.45]),np.array([.6,.8]),1.2*np.eye(2),('seed',0));out=run_dynamic_graph_aims([b],[1+0j],dt=2e-4,steps=20,spawn_threshold=1e-6,overlap_block=.9,max_basis=2,store_every=2)
    assert len(out['events'])==1 and len(out['final_basis'])==2
    norms=np.array([r['norm'] for r in out['records']]);nodes=np.array([r['graph']['nodes'] for r in out['records']]);assert np.all(np.diff(nodes)>=0) and nodes[-1]>nodes[0] and np.max(np.abs(norms-1))<5e-5
    assert abs(out['final_coefficients'][1])>1e-10
