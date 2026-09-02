# v0.5 Build Validation Report

Validated on 2026-08-12.

## Compilation

All Python source files compiled successfully.

## Automated tests

```text
[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m                   [100%][0m
[32m[32m[1m54 passed[0m[32m in 2.23s[0m[0m
```

The v0.5 additions are covered by tests for:

- Cartesian electronic-structure point validation;
- generalized-coordinate gradient/NAC projection;
- projected nuclear mass matrices;
- analytic Gaussian gradient matrix elements;
- analytic Gaussian kinetic matrix elements for a non-diagonal mass matrix;
- Hermiticity of the local adiabatic Gaussian Hamiltonian;
- the moving-basis overlap identity;
- generalized-mass NAC momentum rescaling;
- backend-driven spawning and coupled amplitude transfer;
- deterministic disk caching;
- the explicit PySCF SA-CASSCF call contract;
- PySCF `state=(ket,bra)` to internal `d[i,j]=<i|grad j>` conversion;
- explicit `mult_ediff=False` dynamics calls;
- configured ETF forwarding;
- optional scaled-NAC diagnostic calls;
- explicit RHF/ROHF selection;
- clear behavior when PySCF is absent.

## Representative backend-independent v0.5 outputs

### Coordinate projection

```text
Projected generalized-coordinate point
energies: [0.  0.1]
gradients_q: [ 0.03 -0.02]
d_01/dq: 0.4
mass matrix (electron masses): [[1822.88848621]]
```

### Gridless backend-driven Gaussian dynamics

```text
v0.5 gridless backend-driven Gaussian dynamics
spawn events: [{'step': 1, 'time': 0.0002, 'parent_index': 0, 'new_index': 1, 'target_state': 0}]
basis size: [1 2 2 2 2 2 2 2 2 2 2]
final populations: [3.88661401e-08 9.99999961e-01]
max norm drift: 3.3306690738754696e-16
```

## PySCF runtime status

PySCF is not installed in the build environment used for this release.

Therefore:

- the PySCF backend source was written against the current official PySCF
  SA-CASSCF/gradient/NAC interfaces;
- its call semantics are exercised by a deterministic fake-PySCF integration test;
- the real LiH PySCF examples are included but were not numerically executed here.

Before using the backend for research, install PySCF and run:

```bash
python examples/12_pyscf_lih_sacasscf.py
python examples/14_pyscf_lih_projected_stretch.py
```

Then perform the molecule-specific electronic-structure validation described in
`V05_VALIDATION.md`.
