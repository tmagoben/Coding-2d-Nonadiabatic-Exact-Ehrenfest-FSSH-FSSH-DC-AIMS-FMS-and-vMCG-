# v0.25.1 build validation

Release target: CPython 3.12 on Linux x86-64 with the pinned PySCF 2.13.1 optional
runtime used by inherited gates.

Required release checks:

1. v0.25.1 analytic/TDVP core tests pass.
2. v0.25.1 deterministic evidence passes 55/55 gates.
3. v0.25.1 cumulative benchmark passes 535/535 gates.
4. The full inherited test suite passes in the pinned runtime.
5. Source archive and wheel are built from a cache-free tree.
6. Both artifacts install/import as version 0.25.1 and pass packaged smoke checks.
7. SHA-256 hashes and final file inventory are recorded with the artifacts.

Canonical evidence files:

- `results/v0251_multigaussian_tdvp_evidence.json`
- `results/v0251_multigaussian_tdvp_campaign.json`

Executed full-suite result in the pinned local environment:

- `480 passed in 717.13s (0:11:57)`
- cumulative acceptance: `535/535`
- v0.25.1 core: `19/19`
- v0.25.1 deterministic validation: `55/55` gates

Artifact names, sizes, and hashes are recorded in the adjacent release artifact
inventory after the final cache-free build.
