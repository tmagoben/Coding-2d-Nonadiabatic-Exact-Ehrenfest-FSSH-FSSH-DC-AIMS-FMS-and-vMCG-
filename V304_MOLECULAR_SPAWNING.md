# v0.30.4 — complete-complex spawning scores and molecular checkpoint replay

This development iteration repairs the recovered curved spawning oracle and binds
every generated candidate to a replayable full-Hamiltonian checkpoint. The stable
package remains v0.27.0. No learned model or adaptive trajectory is enabled.

## Corrected candidate space

The v0.30.3 branch still contained the original real-only spawning score. An earlier
local v0.29.6 correction had not reached this branch. The existing TDVP state has
complex electronic coefficients, so each newborn column `K` contributes both `K`
and `iK` as real parameter directions. Omitting the latter loses legitimate defect
reduction and can change the score under a complex electronic gauge transformation.

The repaired v0.29.2 and v0.29.3 evaluators realify the complete candidate space,
project it against the retained current tangent, and compute its singular-value
decomposition. Squared singular values retain the original metric cutoff units.
This also avoids subtracting nearly equal Gram matrices when testing novelty.
The Schur rank now ranges up to `2 * nstate`; it is not an electronic-state count.

The independent reference separately solves the old and augmented real least-squares
systems. It includes both coefficient directions and subtracts the best old-space
residual even when the supplied defect is not perfectly orthogonal to that space.
Additional regression tests cover purely imaginary defects, duplicate spaces,
complex constant and coordinate-dependent gauges, coordinate rotation and packet
permutation. Non-finite snapshots and invalid supplied Wilson transports fail closed.

## Molecular records and replay

`evaluate_molecular_spawn_checkpoint_v304` rebuilds the complete nuclear-kinetic plus
electronic `H psi` and TDVP residual. It records every deterministic position and
momentum candidate, including rejected candidates, with:

- The exact Gaussian state, family/trajectory/step identity and stable packet IDs.
- Hamiltonian and molecular bridge fingerprints, tangent step and oracle settings.
- Candidate geometry, corrected residual gain, retained rank, condition diagnostics
  and all residual-oracle rejection reasons.
- Per-record and complete-checkpoint SHA-256 fingerprints.

Replay reconstructs the state and settings, recomputes the Hamiltonian residual and
every candidate, and compares the entire checkpoint exactly. Re-signing an altered
score or an incomplete candidate list does not make it pass. Exact replay is a
same-runtime contract; physical gauge comparisons use stated numerical tolerances.
Repeated Wilson sections are cached only inside one checkpoint evaluation.

The package-level oracle schema advances from `gnd-spawn-oracle-v1` to
`gnd-spawn-oracle-v2`, with implementation ID
`full-complex-coefficient-tangent-v0.30.4`. The configuration fingerprint also binds
that ID. Old score labels must be recomputed. The independent top-level v0.29.4
`SpawningRecord` schema is unchanged.

## Fresh numerical audit

Run from the repository root with the single-thread environment used in CI:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  python examples/161_recompute_v0304_molecular_spawning.py
```

The audit uses the v0.30.3 analytic singlet–triplet fixture on its 13-point grid.
It evaluates four candidates at both `t=0` and `t=0.05` after one implicit-midpoint
step. Each physical state is transformed into the coordinate-dependent electronic
gauge, then `H psi`, the Jacobian and candidate scores are independently rebuilt.
Using the same physical state isolates the local oracle from trajectory time error.

The score and normalized-gain comparison requires absolute error at most
`3e-10 + 3e-6 * max(abs(a), abs(b))`. Discrete residual decisions, reasons and ranks
must agree exactly. The campaign also checks an independent complete-complex
least-squares reference, nonzero kinetic action, full candidate coverage, exact
moving-frame replay, and a deliberately high-threshold rejection control. Fresh
evidence is written to `results/v0304_molecular_spawning_evidence.json`; CI runs a
separate recomputation and uploads its own result.

The committed local audit passes **13/13 gates** on Python 3.12.14, NumPy 2.5.2
and SciPy 1.18.0. It includes eight physical candidate comparisons (four at each
time), their moving-frame counterparts and four high-threshold rejections.

| Measurement | Maximum observed value |
| --- | ---: |
| Absolute disagreement with independent least squares | 4.34e-19 |
| Absolute score difference between electronic gauges | 1.19e-13 |
| Residual gain missed by the old real-only score | 2.98e-4 |

The focused oracle/molecular bridge/replay tests pass 51 checks, and 16 additional
evidence controls reject corrupted scores, missing candidates, false admission,
failed replay and inconsistent summary columns. The regular wheel builds and its
new modules import outside the repository checkout. Cumulative regression and
fresh hosted campaign results are recorded on the development pull request.

The first hosted cumulative run exposed two pre-existing roundoff-sensitive
assertions: a normalized norm of `1.0000000000000002` exceeded a `2e-16` bound,
and an SVD equation residual of `2.46e-15` exceeded a `2e-15` bound. These tests
now use dimension-scaled float64 roundoff budgets, check that normalization preserves
the electronic direction, and reconstruct the metric equation residual directly. A perturbed
velocity is required to fail the latter bound. The normalization and SVD production
algorithms, null-space checks and scientific campaign gates are unchanged.

## Admission and dataset boundary

These are **residual-screening records**, with `final_admission=null` and
`dataset_certified=false`. Passing the score is not a completed spawn event.
No projection, energy-jump, lifecycle or newborn-activation outcome is invented.
Likewise, a candidate spanning the full electronic spinor has no single target
electronic state. The current top-level `OracleDecision` requires one for `spawn`;
silently setting it to zero would corrupt the scientific meaning of the label.

The remaining adapter must version that distinction explicitly, include final
admission receipts, and audit physical invariant features with numerical tolerances.
In particular, raw parameter-metric condition numbers are not certified as invariant
under coordinate-dependent electronic gauge changes. Gauge variants and all nearby
candidate evaluations must share their trajectory's split. A multi-trajectory,
physically audited campaign is still required before freezing training data.

This audit does not establish continuum convergence, basis completeness, live
ab-initio SOC dynamics, adaptive topology, or a certified training dataset.
