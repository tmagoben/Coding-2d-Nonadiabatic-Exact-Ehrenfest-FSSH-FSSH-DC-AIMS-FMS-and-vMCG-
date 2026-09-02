import numpy as np

from gaussian_dynamics.molecular_backend import MolecularGeometry
from gaussian_dynamics.pyscf_backend_v05 import PySCFSACASSCFConfig
from gaussian_dynamics.pyscf_tracked_backend_v06 import (
    PySCFTrackedSACASSCFBackend,
)

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
    charge=0,
    spin=0,
    scf_reference="RHF",
    warm_start_mo=True,
    use_etfs=False,
    verbose=0,
)

backend=PySCFTrackedSACASSCFBackend(
    config,
    minimum_overlap=0.4,
    minimum_score_margin=0.02,
    ambiguity_policy="raise",
)

distances=np.linspace(2.8,3.2,5)

print("Tracked LiH SA-CASSCF scan")
print("---------------------------")

for R in distances:
    geometry=MolecularGeometry(
        ("Li","H"),
        np.array([
            [0.0,0.0,0.0],
            [0.0,0.0,R],
        ]),
    )

    point=backend.evaluate(geometry)

    print(f"\nR = {R:.6f} bohr")
    print("tracked energies:",point.energies)
    print(
        "tracked->raw:",
        point.metadata.get(
            "permutation_tracked_to_raw",
            list(range(config.nstates)),
        ),
    )
    print(
        "assigned |overlap|:",
        point.metadata.get(
            "assigned_overlap_magnitudes",
            [1.0]*config.nstates,
        ),
    )
    print(
        "overlap unitarity defect:",
        point.metadata.get("state_overlap_unitarity_defect",0.0),
    )
