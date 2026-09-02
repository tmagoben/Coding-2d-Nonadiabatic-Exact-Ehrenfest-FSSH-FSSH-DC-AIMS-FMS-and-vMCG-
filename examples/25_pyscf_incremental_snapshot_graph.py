import numpy as np
from gaussian_dynamics.molecular_backend import MolecularGeometry
from gaussian_dynamics.pyscf_backend_v05 import PySCFSACASSCFConfig
from gaussian_dynamics.pyscf_tracked_backend_v06 import PySCFTrackedSACASSCFBackend
from gaussian_dynamics.incremental_snapshot_graph import IncrementalSnapshotGaugeGraph
try:
    import pyscf
except ImportError:
    raise SystemExit("PySCF is not installed. Install with: pip install -e '.[pyscf]'")
config=PySCFSACASSCFConfig(basis='sto-3g',ncas=2,nelecas=2,nstates=2,weights=(0.5,0.5),warm_start_mo=True,use_etfs=False,verbose=0)
backend=PySCFTrackedSACASSCFBackend(config)
builder=IncrementalSnapshotGaugeGraph(2)
previous=None
for k,R in enumerate([2.9,3.0,3.1]):
    geom=MolecularGeometry(('Li','H'),np.array([[0.,0.,0.],[0.,0.,R]]))
    point,snapshot=backend.evaluate_raw_with_snapshot(geom)
    node=('LiH',k)
    builder.add_cartesian_point(node,snapshot,point,connect_to=[] if previous is None else [previous])
    previous=node
print('nodes:',len(builder.graph.nodes))
print('edges:',len(builder.graph.edges()))
for row in builder.edge_diagnostics():
    print(row)
