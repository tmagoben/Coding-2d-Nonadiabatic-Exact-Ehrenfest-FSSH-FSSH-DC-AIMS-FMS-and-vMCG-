import numpy as np
from gaussian_dynamics.dynamic_graph_aims import DynamicGraphTBF,run_dynamic_graph_aims

parent=DynamicGraphTBF(0,1,np.array([0.55,0.45]),np.array([0.6,0.8]),1.2*np.eye(2),('seed',0))
out=run_dynamic_graph_aims([parent],[1.0+0j],dt=2e-4,steps=50,spawn_threshold=1e-6,overlap_block=0.9,max_basis=2,store_every=5)
print('Spawn events:',out['events'])
print('Final coefficients:',out['final_coefficients'])
print('Final record:',out['records'][-1])
print('Final graph nodes:',len(out['graph'].nodes))
print('Final graph edges:',len(out['graph'].edges()))
