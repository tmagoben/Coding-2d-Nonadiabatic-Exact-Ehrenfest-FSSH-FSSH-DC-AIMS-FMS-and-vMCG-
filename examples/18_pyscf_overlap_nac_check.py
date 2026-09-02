import numpy as np

from gaussian_dynamics.molecular_backend import MolecularGeometry
from gaussian_dynamics.pyscf_backend_v05 import PySCFSACASSCFConfig
from gaussian_dynamics.pyscf_tracked_backend_v06 import (
    PySCFTrackedSACASSCFBackend,
)
from gaussian_dynamics.overlap_transport import directional_nac_from_overlap

try:
    import pyscf
except ImportError:
    raise SystemExit(
        "PySCF is not installed. Install with: pip install -e '.[pyscf]'"
    )

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

backend=PySCFTrackedSACASSCFBackend(
    config,
    minimum_overlap=0.4,
    minimum_score_margin=0.02,
)

R0=3.0
dR=1.0e-3

def geom(R):
    return MolecularGeometry(
        ("Li","H"),
        np.array([
            [0.0,0.0,0.0],
            [0.0,0.0,R],
        ]),
    )

p0=backend.evaluate(geom(R0))
p1=backend.evaluate(geom(R0+dR))

encoded=p1.metadata["state_overlap_matrix_tracked_gauge"]
O=np.array([
    [complex(re,im) for re,im in row]
    for row in encoded
])

d_overlap=directional_nac_from_overlap(O,dR)

# The path moves only the H atom along +z.
d_pyscf=p0.nac_cart[:,:,1,2]

print("Tracked overlap matrix:")
print(O)

print("\nDirectional NAC from overlap / bohr^-1:")
print(d_overlap)

print("\nAnalytic PySCF NAC projected onto H-z / bohr^-1:")
print(d_pyscf)

print("\nDifference:")
print(d_overlap-d_pyscf)

print(
    "\nThis is a finite-step consistency diagnostic. "
    "Decrease dR and verify convergence rather than expecting exact equality "
    "at one finite displacement."
)
