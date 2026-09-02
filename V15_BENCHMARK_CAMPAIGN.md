# v0.15 Cost-Aware Cache Campaign

## 1. Purpose

The campaign asks a narrow question:

> Can the v0.14 adaptive TDSE-defect dynamics be reproduced while eliminating repeated
> Gaussian-pair linear algebra and making the growth decision explicitly aware of
> incremental computational cost?

The release therefore keeps the same analytic LVC physics and the same initial
10-Gaussian residual-selected state.

## 2. Initial basis

The initial basis is built exactly as in v0.14:

```text
64 x 64 initial projection grid
726-member deterministic Gaussian dictionary
top-30 residual shortlist
initial reduced-density screening
final initial basis size = 10
```

## 3. Adaptive settings

```text
dt = 0.005
steps = 120
final time = 0.6

defect grid = 40 x 40
defect interval = 10 steps

add threshold = 0.020
remove threshold = 0.006

minimum capture fraction = 0.003
minimum cost-aware utility = 0.15

cost horizon = 10 steps
conditioning penalty weight = 0.15
residual shortlist = 8

minimum basis = 8
maximum basis = 11
```

## 4. Matrix architecture

At each endpoint:

```text
construct one GaussianPairCache
        ↓
one canonical pair solve per i<=j
        ↓
assemble S and H from cached moments
```

At each midpoint:

```text
construct one midpoint GaussianPairCache
        ↓
one canonical pair solve per i<=j
        ↓
assemble ordered T_ij from cached/reversed pair moments
```

At a defect checkpoint:

```text
reuse endpoint GaussianPairCache
        ↓
build endpoint T without new pair factorizations
        ↓
solve projected Cdot
        ↓
reconstruct Psi and Psidot
        ↓
apply independent FFT-grid H
```

## 5. Candidate workflow

```text
generate energy-conserving local candidates
        ↓
vectorized K x G TDSE-defect ranking
        ↓
residual shortlist
        ↓
exact conditioning with temporary expanded pair caches
        ↓
cost-aware utility
        ↓
accept/reject
```

The accepted candidate cache is reused directly by the incremental matrix expansion.

## 6. Cost model

For each residual-qualified candidate:

```text
benefit = predicted TDSE-defect capture fraction

cost =
    relative endpoint+midpoint pair growth
  + relative Cayley cubic growth
  + 0.25 x relative defect-solve cubic growth
  + conditioning multiplier
```

The selector uses

$$
U_c=f_c/C_c.
$$

No future exact observable is used.

## 7. Physics comparison

The exact projected-state reference and exact target are recomputed as in previous
releases.

The v0.15 result is also compared against the saved v0.14 reference metrics.

This is an invariance test, not only an accuracy test.

## 8. Complexity comparison

The campaign records both:

```text
actual v0.15 pair factorizations
v0.14 factorization-equivalent baseline
```

plus cache requests/hits/reverse views, candidate-search pair work, cubic solve units,
and category timings.

## 9. Reproduction

Run:

```bash
python examples/64_recompute_v015_campaign.py
```

Canonical output:

```text
results/v015_cost_aware_cache_campaign.json
```
