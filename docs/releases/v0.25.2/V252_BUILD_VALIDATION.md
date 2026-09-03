# v0.25.2 build validation

Release target: CPython 3.12 on Linux x86-64 with pinned PySCF 2.13.1 for inherited
runtime gates.

Required checks:

1. Adaptive analytic/TDVP core passes 20 tests.
2. Deterministic adaptive evidence passes 70/70 gates.
3. Cumulative release benchmark passes 630/630 gates.
4. The complete inherited/new pytest suite passes in the pinned runtime.
5. The source archive and wheel are cache-free and integrity tested.
6. An isolated wheel install and extracted-source import report version 0.25.2.
7. Artifact SHA-256 hashes are recorded alongside the deliverables.

Canonical evidence paths:

- `results/v0252_adaptive_multigaussian_tdvp_evidence.json`
- `results/v0252_adaptive_multigaussian_tdvp_campaign.json`

Full-suite receipt:

- command: `python -m pytest -q -p no:cacheprovider`
- deterministic runtime controls: one thread for OpenBLAS, OMP, MKL, BLIS,
  NumExpr, and vecLib
- result: **507 passed**
- wall time reported by pytest: **1000.34 s (16 min 40 s)**
- failures, errors, skips, and xfails: **0**

Artifact verification is recorded alongside the final SHA-256 manifest after the
release build.
