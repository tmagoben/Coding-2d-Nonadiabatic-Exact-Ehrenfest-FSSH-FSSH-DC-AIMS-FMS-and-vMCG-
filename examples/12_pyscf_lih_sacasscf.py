import numpy as np

from gaussian_dynamics.molecular_backend import MolecularGeometry
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

geometry=MolecularGeometry(
    symbols=("Li","H"),
    coords_bohr=np.array([
        [0.0,0.0,0.0],
        [0.0,0.0,3.0],
    ]),
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
    scf_conv_tol=1e-10,
    mc_conv_tol=1e-9,
    mc_conv_tol_grad=1e-5,
    use_etfs=False,
    compute_scaled_nac=True,
    warm_start_mo=False,
    verbose=0,
)

backend=PySCFSACASSCFBackend(config)
point=backend.evaluate(geometry)

print("PySCF version:",point.metadata["pyscf_version"])
print("SA-CASSCF energies / Eh:")
print(point.energies)

print("\nState gradients / Eh bohr^-1:")
print(point.gradients_cart)

print("\nInternal d_01 = <0|grad_R 1> / bohr^-1:")
print(point.nac_cart[0,1])

print("\nNAC antisymmetry residual:")
print(np.linalg.norm(point.nac_cart[0,1]+point.nac_cart[1,0]))

print("\nPySCF mult_ediff diagnostic for pair 0,1:")
print(point.scaled_nac_cart[0,1])

print("\nMetadata:")
for key in [
    "backend",
    "basis",
    "ncas",
    "nelecas",
    "state_average_weights",
    "use_etfs",
    "nac_internal_convention",
    "pyscf_request_for_internal_dij",
]:
    print(f"{key}: {point.metadata[key]}")
