# v0.27.0 Build Validation

Release date: 2026-08-25

## Required checks

1. compile every Python source file;
2. run the v0.27.0 correlated core and lifecycle tests;
3. recompute the 100-gate scientific evidence;
4. recompute the cumulative 960-gate campaign, including all 825 inherited gates;
5. run the complete test suite with the pinned PySCF 2.13.1 environment;
6. build source and wheel artifacts;
7. install the wheel into an isolated directory and rerun package smoke tests;
8. verify deterministic archive manifests and SHA-256 checksums.

## Frozen environment

- Python: 3.12
- NumPy: 2.5.2 in the pinned PySCF validation environment
- SciPy: 1.18.0 in the pinned PySCF validation environment
- PySCF optional extra: 2.13.1

The core correlated release does not require PySCF at runtime. PySCF is retained to
exercise all inherited molecular-backend gates.

## Scientific acceptance

- Correlated TDVP and lifecycle core tests: **18 passed in 3.14 s**.
- Deterministic correlated-width evidence: **100/100 gates passed**.
- Evidence fingerprint:
  `d3a7509c517b2c44cb5e058d80fabf29d549142589c7b3875df17372611fc323`.
- Cumulative release campaign: **960/960 gates passed**.
- Inherited gates: 825; new scientific gates: 100; new adversarial/core gates: 35.
- Cumulative campaign JSON SHA-256:
  `6aae2619cc4e2f215b138a38c10ce7566bc43304c4e1256a1e496f0f1c5c3800`.

The campaign and complete suite were launched with OpenBLAS, OMP, MKL, BLIS,
NumExpr, and vecLib fixed to one thread before importing NumPy or PySCF. The
release-locked runtime identity passed exactly: CPython 3.12.13, NumPy 2.5.2,
SciPy 1.18.0, h5py 3.16.0, and PySCF 2.13.1.

## Complete regression receipt

- command: `python -m pytest -q -p no:cacheprovider`
- result: **585 passed**
- wall time reported by pytest: **1378.80 s (22 min 58.80 s)**
- failures, errors, skips, and xfails: **0**
- Python sources compiled independently in memory: **538**

## Artifact receipts

- wheel: `gaussian_nadyn-0.27.0-py3-none-any.whl`
- wheel SHA-256:
  `c63b8bafa9e58823480426008fb907631c05a97a94590b67a0a8aa130939af1d`
- two clean offline wheel builds with `SOURCE_DATE_EPOCH=1787616000` were
  byte-for-byte identical
- wheel inventory: 219 entries; all four v0.27.0 modules are present; no bytecode
  or cache entries are present
- isolated wheel installation loaded v0.27.0 outside the source tree; the
  correlated metric had rank 14 and nullity 0, and an implicit-midpoint step had
  zero reported norm change, `2.17e-19` energy change, and a positive width matrix
- deterministic source ZIP inventory: 1005 files under one fixed root, with fixed
  timestamps and no build, egg-info, pytest, or bytecode caches
- two source ZIP builds were byte-for-byte identical
- extracted-source import, campaign inspection, and 23 focused model, validation,
  and packaging tests passed

The adjacent `SHA256SUMS.txt` is the authoritative digest manifest for the final
wheel and source ZIP.
