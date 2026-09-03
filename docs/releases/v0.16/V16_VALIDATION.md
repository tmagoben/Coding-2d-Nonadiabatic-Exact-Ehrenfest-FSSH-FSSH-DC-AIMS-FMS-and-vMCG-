# v0.16 Validation Contract

v0.16 introduces an actual sparse approximation, so validation separates:

1. exact algebraic properties of the locality machinery;
2. approximation error introduced by dropped edges;
3. physical benchmark error;
4. scaling behavior.

## Exact automated checks

The cumulative test suite verifies:

- the conservative position-overlap bound is never below the exact overlap on
  randomized SPD Gaussian pairs;
- far pairs can be rejected without an exact pair solve;
- enter/exit hysteresis retains and removes edges correctly;
- full-graph sparse S/H equals dense v0.15 S/H;
- full-graph sparse T and Cayley propagation equal dense v0.15 results;
- the geometry-cache electronic cost model distinguishes hits from new points;
- local-degree cost uses predicted degree rather than total basis size;
- provider cost can change candidate ordering;
- short sparse adaptive propagation performs a real enrichment and preserves norm;
- release acceptance rejects an effectively dense graph when sparsity is required.

## Release thresholds

```text
initial reduced-density error       <= 0.035
projected-state dynamics error      <= 0.001
target density error                <= 0.035
population L2 error                 <= 0.03
coherence phase error               <= 0.0035 rad
generalized norm drift              <= 1e-4
condition number                    <= 5e3

average compact-benchmark sparsity  >= 4%
pair reduction vs v0.15             >= 4%
final rho difference vs v0.15       <= 0.0015

N=80 chain pair reduction            >= 90%
N=80 chain edge fraction             <= 0.08

fitted local-edge exponent           <= 1.20
fitted dense-pair exponent           >= 1.80

final sparse S relative error        <= 0.01
final sparse H relative error        <= 0.01
```

## Measured release values

```text
projected-state dynamics error:
0.00013361460054812487

target density error:
0.03333954068459046

population error:
0.028199413658981914

coherence phase error:
0.0029095064228781115

norm drift:
2.0053154308197207e-06

compact graph sparsity:
0.053816793893129766

pair reduction vs v0.15:
0.04478468899521526

final rho difference vs v0.15:
0.0001631265880682423

final sparse S relative error:
0.005191742661052565

final sparse H relative error:
0.003962632349871911

N=80 pair reduction:
0.9268518518518518

N=80 edge fraction:
0.04968354430379747

active-edge exponent:
1.0425836916313385

dense-pair exponent:
1.9737662900529327
```

All configured release checks pass.

## Sparse operator audit

The overlap cutoff rigorously bounds the dropped overlap magnitude, but does not by
itself provide a universal Hamiltonian error bound.

Therefore the release endpoint is rebuilt densely and compared with the sparse matrix.

The dense audit is intentionally expensive and is not performed on every time step.

## Physical comparison against v0.15

The final reduced density differs from the dense v0.15 reference by

$$
0.0001631265880682423
$$

in Frobenius norm.

This is below the release threshold while still being large enough to show that v0.16
is a genuine sparse approximation rather than merely a different storage format.

## Timing policy

Wall-clock speedups are diagnostic only.

Release acceptance uses operation counts, graph sparsity, dense-reference agreement,
and physical observables.

## PySCF scope

The electronic cost model is validated as an abstract cache-cost interface.

No real PySCF wall-time predictor is claimed in v0.16.
