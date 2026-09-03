# Optional PySCF SA-CASSCF Provider

The optional adapter follows the public PySCF interfaces for:

- state-averaged CASSCF;
- state-specific SA-CASSCF nuclear gradients;
- SA-CASSCF analytical nonadiabatic couplings.

## Minimal conceptual use

```python
from gaussian_dynamics.pyscf_provider import PySCFStateAveragedCASSCFProvider

def h2_coordinate(R):
    atoms = [
        ("H", (0.0, 0.0, -R/2)),
        ("H", (0.0, 0.0,  R/2)),
    ]
    tangent = [
        [0.0, 0.0, -0.5],
        [0.0, 0.0,  0.5],
    ]
    return atoms, tangent

provider = PySCFStateAveragedCASSCFProvider(
    geometry_builder=h2_coordinate,
    basis="sto-3g",
    ncas=2,
    nelecas=2,
    nstates=2,
)

point = provider.evaluate(1.5)  # q in bohr
```

This example is a software-interface demonstration, not a recommended production
electronic-structure level for photochemical dynamics.

## Important convention

PySCF documents its SA-CASSCF NAC state tuple as `(ket, bra)` and returns

`<bra | d(ket)/dR>`.

The adapter stores instead

`nac_q[i,j] = <i | d/dq j>`,

so v0.3 requested PySCF `state=(j,i)` and then projected the Cartesian NAC onto
the generalized coordinate tangent. **v0.23.2 erratum:** real PySCF 2.13.1
overlap finite differences show that the production mapping for this internal
convention is `state=(i,j)`, `mult_ediff=False`, `use_etfs=False`. See
`../v0.23.2/V232_NAC_CONVENTION_ERRATUM.md`.

## Scaled NAC warning

PySCF exposes a `mult_ediff` option that returns the NAC multiplied by the state
energy difference. This v0.3 contract does **not** use that option because downstream
dynamics expects the actual derivative coupling.

## Runtime-validation status

The core v0.3 provider architecture is tested with analytic and tabulated providers.
The PySCF adapter is optional and is not required for the core test suite. Before using
it for research calculations, run a backend-specific validation set for the chosen
molecule, active space, state manifold, and PySCF version.
