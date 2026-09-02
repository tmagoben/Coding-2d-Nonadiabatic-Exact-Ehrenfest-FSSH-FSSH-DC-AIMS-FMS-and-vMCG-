# v0.11 Validation Contract

Version 0.11 is accepted only if all earlier regression tests remain passing and the
new basis-completeness mathematics passes independent checks.

## Unequal-width Gaussian algebra

The test suite verifies:

- exact overlap against two-dimensional numerical quadrature;
- exact gradient matrix element against quadrature;
- exact kinetic matrix element against quadrature with a non-diagonal mass matrix;
- reduction to the old equal-width complex centroid;
- reduction to the arithmetic saddle point for equal widths;
- moving-basis derivative including $\dot A$ against finite differences.

## Spawn-candidate search

For every accepted synthetic candidate:

- parent and child classical energies agree within tolerance;
- the search is deterministic;
- candidates are ranked by decreasing score;
- sibling selection does not return identical phase-space/width copies;
- redundant target-state candidates are blocked.

## v0.11 managed dynamics

A short strong-coupling regression must show:

- at least one optimized spawning event;
- a finite positive-definite reduced electronic density;
- normalized reduced density trace;
- controlled generalized norm;
- lineage metadata for every created TBF.

## Basis completeness diagnostics

The tests verify:

- canonical coefficient norm equals $C^\dagger SC$;
- participation ratio lies within the retained basis dimension;
- generation histograms are deterministic;
- width-diversity diagnostics use determinant ratios consistently.

## Full cumulative regression

All v0.1-v0.10 tests remain required.

A v0.11 improvement is not allowed to break:

- exact FFT references;
- PySCF API contracts;
- many-electron state tracking;
- gauge graph covariance;
- Wilson-loop topology;
- moving-basis metric compatibility;
- SPA0/SPA1 tests;
- v0.10 benchmark/observable tests.

## Strong-CI release benchmark

The release campaign compares:

1. exact TDSE;
2. v0.10 baseline;
3. v0.11 basis ladder $N_{\max}=2,4,6,8,10$;
4. no-position-optimization ablation;
5. fixed-width ablation;
6. single-child-per-event ablation.

The following observables are reported separately:

- global-diabatic population error;
- full reduced-density Frobenius error;
- purity error;
- generalized norm drift;
- maximum overlap condition number.

The default acceptance thresholds intentionally distinguish diagonal population
accuracy from full-density/coherence accuracy.

A run can therefore pass the population criterion while failing the full-density
criterion.

That outcome must be reported as partial convergence, not full validation.
