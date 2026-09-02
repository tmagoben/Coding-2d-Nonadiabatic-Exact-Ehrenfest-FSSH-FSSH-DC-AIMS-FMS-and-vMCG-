import numpy as np

from gaussian_dynamics.molecular_backend import (
    MolecularGeometry,
    LinearGeometryMap,
    GeneralizedCoordinateProvider,
)
from gaussian_dynamics.pyscf_backend_v05 import (
    PySCFSACASSCFConfig,
    PySCFSACASSCFBackend,
)

try:
    import pyscf
except ImportError:
    raise SystemExit(
        "PySCF is not installed. Install with: pip install -e '.[pyscf]'"
    )

R0=3.0

reference=MolecularGeometry(
    ("Li","H"),
    np.array([
        [0.0,0.0,0.0],
        [0.0,0.0,R0],
    ]),
)

backend=PySCFSACASSCFBackend(
    PySCFSACASSCFConfig(
        basis="sto-3g",
        ncas=2,
        nelecas=2,
        nstates=2,
        weights=(0.5,0.5),
        warm_start_mo=True,
        use_etfs=False,
        verbose=0,
    )
)

# First point supplies the same PySCF masses used by the dynamics projection.
raw=backend.evaluate(reference)
m_li,m_h=raw.masses_amu
mtot=m_li+m_h

# COM-preserving bond-stretch tangent dR_cart/dq.
modes=np.zeros((1,2,3))
modes[0,0,2]=-m_h/mtot
modes[0,1,2]= m_li/mtot

geomap=LinearGeometryMap(
    reference.symbols,
    reference.coords_bohr,
    modes,
)
provider=GeneralizedCoordinateProvider(backend,geomap)

point=provider.evaluate(np.array([0.0]))

print("LiH projected bond-stretch SA-CASSCF point")
print("PySCF version:",point.metadata["pyscf_version"])
print("energies / Eh:",point.energies)
print("dE/dq / Eh bohr^-1:",point.gradients_q[:,0])
print("d_01/dq / bohr^-1:",point.nac_q[0,1,0])
print("effective stretch mass / electron masses:",point.mass_matrix_q_au[0,0])
