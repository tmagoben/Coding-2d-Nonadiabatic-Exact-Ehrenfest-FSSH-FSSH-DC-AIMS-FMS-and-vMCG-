# Explicit PySCF SA-CASSCF Backend

This document describes the **actual backend implementation** in

```text
gaussian_dynamics/pyscf_backend_v05.py
```

rather than a pseudocode interface.

---

## 1. Installation

The core repository does not require PySCF.

To enable the backend:

```bash
pip install -e ".[pyscf]"
```

or install PySCF in the active environment and then install the project normally.

---

## 2. Minimal backend construction

```python
from gaussian_dynamics.pyscf_backend_v05 import (
    PySCFSACASSCFConfig,
    PySCFSACASSCFBackend,
)

config = PySCFSACASSCFConfig(
    basis="cc-pvdz",
    ncas=4,
    nelecas=4,
    nstates=2,
    weights=(0.5, 0.5),

    charge=0,
    spin=0,
    scf_reference="RHF",

    scf_conv_tol=1e-10,
    mc_conv_tol=1e-9,
    mc_conv_tol_grad=1e-5,

    use_etfs=False,
    compute_scaled_nac=False,
    warm_start_mo=True,
)

backend = PySCFSACASSCFBackend(config)
```

The active space shown here is an **example**, not a universal recommendation.
Active-space selection is a chemical/scientific convergence problem.

---

## 3. Geometry input

All backend coordinates are explicit **Bohr**:

```python
import numpy as np

from gaussian_dynamics.molecular_backend import MolecularGeometry

geometry = MolecularGeometry(
    symbols=("Li", "H"),
    coords_bohr=np.array([
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 3.0],
    ]),
)
```

Evaluation is then

```python
point = backend.evaluate(geometry)
```

The point contains

```python
point.energies
point.gradients_cart
point.nac_cart
point.masses_amu
point.metadata
```

with shapes

```text
energies        (nstate,)
gradients_cart  (nstate, natom, 3)
nac_cart        (nstate, nstate, natom, 3)
masses_amu      (natom,)
```

---

## 4. What the backend runs

The implementation executes the equivalent of:

```python
mol = gto.M(
    atom=atoms,
    basis=basis,
    charge=charge,
    spin=spin,
    unit="Bohr",
)

mf = scf.RHF(mol)
mf.conv_tol = scf_conv_tol
mf.kernel()

mc = mcscf.CASSCF(mf, ncas, nelecas)
mc = mc.state_average_(weights)

mc.conv_tol = mc_conv_tol
mc.conv_tol_grad = mc_conv_tol_grad
mc.kernel()
```

Then each state gradient is obtained by

```python
mc.nuc_grad_method(state=I).kernel()
```

and each derivative coupling uses the SA-CASSCF NAC interface.

---

## 5. NAC index conversion

The repository convention is

$$
d_{IJ}
=
\langle\Phi_I|\nabla\Phi_J\rangle.
$$

PySCF documents

```text
state=(ket, bra)
```

and returns

$$
\langle\mathrm{bra}|\nabla\,\mathrm{ket}\rangle.
$$

This was the literal interpretation used through v0.23.1:

```python
nac_method.kernel(
    state=(j, i),
    use_etfs=...,
    mult_ediff=False,
)
```

for the internal array

```python
nac_cart[i, j]
```

**v0.23.2 erratum:** direct PySCF 2.13.1 execution and phase-aligned
many-electron overlap derivatives establish `state=(i,j)` as the production
mapping for `nac_cart[i,j]`. The centralized v0.23.2 helper enforces the corrected
mapping with `mult_ediff=False` and `use_etfs=False`. See
`../v0.23.2/V232_NAC_CONVENTION_ERRATUM.md`.

and then sets

```python
nac_cart[j, i] = -nac_cart[i, j]
```

for the real-state convention.

There is no implicit index guessing downstream.

---

## 6. `mult_ediff`

The dynamics always uses

```python
mult_ediff=False
```

because it needs the actual derivative coupling.

If

```python
compute_scaled_nac=True
```

the backend makes a second diagnostic call with

```python
mult_ediff=True
```

and stores that PySCF-specific scaled object in

```python
point.scaled_nac_cart
```

It is **not substituted into the dynamics NAC field**.

---

## 7. ETF choice

The PySCF NAC implementation exposes

```python
use_etfs=True/False
```

for electron translation factors.

v0.5 forwards the configured value explicitly and records it in metadata.

Do not compare NAC data from calculations using different ETF conventions without
noting the distinction.

---

## 8. Warm-start orbitals

If

```python
warm_start_mo=True
```

the converged CASSCF orbital matrix from the preceding call is passed as the initial
orbital guess to the next SA-CASSCF calculation.

This is intended for nearby direct-dynamics geometries.

It does not replace electronic-state tracking.

---

## 9. LiH backend smoke test

`examples/12_pyscf_lih_sacasscf.py` provides a deliberately small backend exercise
using two-state SA-CASSCF(2,2).

Its purpose is to verify:

- PySCF installation;
- SCF convergence;
- SA-CASSCF convergence;
- state-energy extraction;
- state-specific gradients;
- pair NAC extraction;
- the repository's NAC index convention.

It is **not** offered as a quantitatively converged photochemical model.

---

## 10. Projection onto a bond-stretch coordinate

For a one-dimensional molecular coordinate $q$,

$$
R(q)=R_0+Jq.
$$

Build

```python
LinearGeometryMap(...)
```

and wrap the backend:

```python
provider = GeneralizedCoordinateProvider(
    backend,
    geometry_map,
)
```

Then

```python
point_q = provider.evaluate([0.0])
```

returns

$$
E_I(q),
\qquad
\frac{dE_I}{dq},
\qquad
d_{IJ}^{(q)},
\qquad
M_q.
$$

These are the exact objects consumed by the v0.5 Gaussian direct-dynamics layer.

---

## 11. What must still be checked for a real molecule

Do not treat a successful PySCF call as a complete electronic-structure validation.

Before using the backend for scientific dynamics, verify:

1. active orbitals;
2. active electron count;
3. number and character of averaged states;
4. state-average weights;
5. state ordering along the relevant region;
6. orbital/root continuity;
7. gradient finite differences;
8. NAC convention and continuity;
9. basis-set dependence;
10. whether the selected electronic-structure level is adequate for the target
    crossing/seam.

v0.5 makes the backend explicit; it does not automate those chemical decisions.


## 12. End-to-end projected LiH stretch example

`examples/14_pyscf_lih_projected_stretch.py` goes one step beyond the raw backend.

It:

1. evaluates LiH once with the real PySCF backend;
2. reads the PySCF atomic masses;
3. constructs a center-of-mass-preserving Li-H stretch tangent;
4. builds a `LinearGeometryMap`;
5. projects Cartesian state gradients and NAC vectors to the bond coordinate;
6. constructs the corresponding generalized-coordinate nuclear mass.

Thus the printed

```text
energies
dE/dq
d_01/dq
M_q
```

are exactly the quantities expected by the backend-driven Gaussian dynamics.
