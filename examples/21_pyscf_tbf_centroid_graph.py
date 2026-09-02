"""Skeleton for a real PySCF branched TBF/centroid gauge graph.

This intentionally stops after graph construction/diagnostics so users inspect the
electronic overlap network before assembling dynamics matrices.
"""

import numpy as np

from gaussian_dynamics.molecular_backend import MolecularGeometry
from gaussian_dynamics.pyscf_backend_v05 import PySCFSACASSCFConfig
from gaussian_dynamics.pyscf_tracked_backend_v06 import PySCFTrackedSACASSCFBackend
from gaussian_dynamics.pyscf_gauge_graph import (
    build_snapshot_gauge_graph,
    edge_overlap_diagnostics,
    tbf_centroid_edge_pairs,
)

try:
    import pyscf
except ImportError:
    raise SystemExit("PySCF is not installed. Install with: pip install -e '.[pyscf]'")

config=PySCFSACASSCFConfig(
    basis="sto-3g",
    ncas=2,
    nelecas=2,
    nstates=2,
    weights=(0.5,0.5),
    warm_start_mo=True,
    use_etfs=False,
    verbose=0,
)

# Use separate sequential tracker instances for the center calculations in a real
# workflow if they belong to different trajectory histories.  Here we only need raw
# converged snapshots at nearby demonstration geometries.
backend=PySCFTrackedSACASSCFBackend(config,ambiguity_policy="warn")

R_values={
    "t0":2.90,
    "t1":3.00,
    "t2":3.10,
    "c01":2.95,
    "c12":3.05,
    "c02":3.00,
}

snapshots={}
for node,R in R_values.items():
    geometry=MolecularGeometry(
        ("Li","H"),
        np.array([[0.0,0.0,0.0],[0.0,0.0,R]])
    )
    backend.reset_tracking(reset_orbitals=False)
    backend.evaluate(geometry)
    snapshots[node]=backend.previous_snapshot

edge_pairs=tbf_centroid_edge_pairs(
    ["t0","t1","t2"],
    {(0,1):"c01",(1,2):"c12",(0,2):"c02"},
)

graph=build_snapshot_gauge_graph(snapshots,edge_pairs)

print("PySCF TBF-centroid gauge graph")
print("--------------------------------")
for item in edge_overlap_diagnostics(graph):
    print(item)

for cycle in graph.fundamental_cycles("t0"):
    W=graph.wilson_loop(cycle)
    print("\ncycle:",cycle)
    print("Wilson eigenvalues:",np.linalg.eigvals(W))
