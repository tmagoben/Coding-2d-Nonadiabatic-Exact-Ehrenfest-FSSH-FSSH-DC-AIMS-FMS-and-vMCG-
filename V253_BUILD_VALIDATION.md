# v0.25.3 build validation

Release target: CPython 3.12 on Linux x86-64 with pinned PySCF 2.13.1 for inherited
runtime gates.

Required checks:

1. Controlled-basis core passes 21 tests.
2. Deterministic lifecycle evidence passes 60/60 gates.
3. Cumulative release benchmark passes 715/715 gates.
4. The complete inherited/new pytest suite passes in the pinned runtime.
5. Source archive and wheel are cache-free and integrity tested.
6. An isolated wheel install and extracted-source import report version 0.25.3.
7. Both release artifacts rebuild byte-for-byte and have recorded SHA-256 hashes.

Full-suite receipt:

- command: `python -m pytest -q -p no:cacheprovider`
- deterministic runtime controls: one thread for OpenBLAS, OMP, MKL, BLIS,
  NumExpr, and vecLib
- result: **535 passed**
- wall time reported by pytest: **1010.90 s (16 min 50 s)**
- failures, errors, skips, and xfails: **0**

Artifact integrity and reproducibility are recorded in the SHA-256 manifest after
sealing.
