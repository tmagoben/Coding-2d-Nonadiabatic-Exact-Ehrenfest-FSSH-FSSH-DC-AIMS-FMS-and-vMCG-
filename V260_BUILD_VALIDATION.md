# v0.26.0 build validation

Release target: CPython 3.12 on Linux x86-64 with the hash-locked PySCF 2.13.1
runtime recorded in `requirements-pyscf-v260-linux-x86_64-py312.txt`.

## Scientific acceptance

- Focused multidimensional model, exact-grid, TDVP, lifecycle, and evidence tests:
  **25 passed**.
- Deterministic multidimensional evidence: **80/80 gates passed**.
- Cumulative release campaign: **825/825 gates passed**.
- Inherited gates: 715; new scientific gates: 80; new adversarial/core gates: 30.
- Evidence fingerprint:
  `1d2a0e483bb0abc533573acd397d8c8061cd43ba7a59fff04616c102e1a1fe50`.
- Cumulative campaign JSON SHA-256:
  `b4cb8ee7a5aa38ff7c7fecce4c889c19167f96610175b0cf7678e22b6ac3f24f`.

The campaign was run with OpenBLAS, OMP, MKL, BLIS, NumExpr, and vecLib fixed to
one thread before importing NumPy or PySCF.

## Release sealing checklist

Full-suite receipt:

- command: `python -m pytest -q -p no:cacheprovider`
- deterministic runtime controls: one thread for OpenBLAS, OMP, MKL, BLIS,
  NumExpr, and vecLib
- result: **562 passed**
- wall time reported by pytest: **1181.58 s (19 min 41 s)**
- failures, errors, skips, and xfails: **0**

Artifact receipts:

- wheel: `gaussian_nadyn-0.26.0-py3-none-any.whl`
- wheel SHA-256:
  `ebf1916de5d5949e2af4622b40eadbb5a099d8d8d9fedf8a1b8e89ca8d6a1e87`
- two clean wheel builds with `SOURCE_DATE_EPOCH=1787529600` were byte-for-byte
  identical
- wheel inventory: 215 entries, all five v0.26.0 modules present, no bytecode or
  cache entries
- fresh isolated wheel import: version 0.26.0 loaded from the new virtual
  environment; a 2D implicit-midpoint TDVP smoke step completed
- the deterministic source ZIP rebuild is byte-for-byte identical, excludes build,
  egg-info, pytest, and bytecode caches, and imports as version 0.26.0 after extraction

The adjacent `SHA256SUMS.txt` is the authoritative digest manifest for both final
release artifacts.
