# v0.24.0 build validation

Release date: 2026-08-21

## Canonical environment

- CPython 3.12.13 on Linux x86-64, little endian
- NumPy 2.5.2
- SciPy 1.18.0
- h5py 3.16.0
- PySCF 2.13.1 from the preserved hash-locked wheel runtime
- BLAS/OpenMP-related thread controls fixed to one
- OpenMolcas: not installed and not executed

## Recorded validation

- Focused v0.24.0 release/packaging tests: **2 passed in 41.03 s**.
- Canonical acceptance campaign: **256/256 gates passed** (208 inherited, 48 new).
- Full cumulative source-tree suite: **404 passed in 346.10 s**.
- Two independent campaign writes were byte-identical; canonical campaign SHA-256:
  `3668e40280e6265e00f058d37944b246e14c4efe7ffe4b9afa8986e93000aa40`.
- Wheel/isolated-install check: not executed. The environment approval-usage quota
  rejected the packaging command before it ran; this is not recorded as a pass.

## Claim boundary

The build validates the protocol, strict parser, artifact trust chain, derivative and
independent-validation schemas, admission logic, and fixture-bound dynamics rehearsal.
The fixture is synthetic and cannot be relabeled or admitted. External/live
molecular-SOC and ab-initio SOC claims remain false. Native OpenMolcas numerical
cross-parsing remains explicitly unimplemented and is a closed admission prerequisite.
