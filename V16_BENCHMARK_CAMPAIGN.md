# v0.16 Sparse Locality Benchmark Campaign

## Purpose

The campaign tests whether an error-controlled local graph can reduce Gaussian pair
work while preserving the v0.15 physical benchmark within a documented sparse
approximation tolerance.

It also separates two regimes:

```text
compact CI packet basis:
    highly overlapping
    only mildly sparse

bounded-locality Gaussian chain:
    finite local degree
    strongly sparse
```

## Primary CI benchmark

```text
dt = 0.005
steps = 120
final time = 0.6
defect grid = 40 x 40

locality enter threshold = 0.03
locality exit threshold = 0.015

add threshold = 0.020
remove threshold = 0.006
minimum capture fraction = 0.003
minimum local utility = 0.08

initial basis = 10
maximum basis = 11
```

The initial projected state is the same residual-selected 10-Gaussian state used by
v0.15.

## Locality workflow

At every endpoint and midpoint:

```text
precompute minimum width eigenvalues
        ↓
safe KD-tree global-radius query
        ↓
pair-specific conservative overlap upper bound
        ↓
exact overlap only for surviving candidates
        ↓
enter/exit hysteresis
        ↓
sparse exact pair matrices on active edges
```

## Sparse TDSE propagation

The projected equation is propagated with sparse CSR/CSC matrices and
`scipy.sparse.linalg.spsolve`.

The physical-grid TDSE defect remains independent of the sparse matrix storage.

## Adaptive growth

The residual shortlist is generated exactly as before.

The cost reranking now depends on:

```text
predicted local graph degree
estimated sparse block growth
conditioning
electronic-structure cache cost
```

rather than only global basis size.

## Dense endpoint audit

After the release propagation, the final sparse S/H matrices are compared with a full
dense pair build.

This audit is part of release validation because an overlap threshold does not
constitute a universal H-matrix error theorem.

## Bounded-locality scaling benchmark

Synthetic chains are evaluated at

```text
N = 20, 40, 80
```

with the same locality thresholds.

For each size the campaign records:

```text
active edges
KD-tree spatial candidates
global screens
exact overlap checks
pair factorizations
dense canonical pairs
sparse H density
dense assembly wall time
sparse assembly wall time
```

A log-log exponent is fit only as a descriptive scaling diagnostic for this chain.

## Electronic-cost demonstration

A geometry-cache model assigns one candidate a cache-hit cost and another a new-point
cost.

The demonstration verifies that provider cost enters the local utility without
changing the analytic LVC release physics.

## Reproduction

Run:

```bash
python examples/70_recompute_v016_campaign.py
```

Canonical output:

```text
results/v016_sparse_locality_campaign.json
```
