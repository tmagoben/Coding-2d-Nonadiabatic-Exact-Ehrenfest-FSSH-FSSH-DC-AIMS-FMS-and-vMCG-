# v0.15 Validation Contract

v0.15 is a performance/architecture release, so it has a stricter invariance
requirement than a release that deliberately changes the physical approximation.

## 1. Cumulative regression

Every retained v0.1-v0.14 test must remain passing.

New v0.15 tests independently check:

```text
cached S/H against the original exact S/H builder
cached T against the original moving-basis T builder
reverse pair identities
incremental add against a full rebuild
incremental subset/prune against a full rebuild
factorization bookkeeping
cost-model algebra
cost-aware conditioning penalty
cost utility gate
v0.15 fixed-basis propagation against v0.14
short real cost-aware enrichment
release acceptance logic
```

## 2. Pair-cache equality

For mixed-width/mixed-momentum test Gaussians, the cached implementation must reproduce
the old analytic matrices:

$$
S_{cache}=S_{reference},
$$

$$
H_{cache}=H_{reference},
$$

$$
T_{cache}=T_{reference}
$$

to tight floating-point tolerance.

## 3. Reverse pair

After computing $(i,j)$, requesting $(j,i)$ must:

```text
not increase canonical pair solve count
conjugate the overlap
conjugate the cross centroid
reuse the covariance
```

## 4. Incremental add

Starting from an $N$-Gaussian matrix, append one Gaussian.

The incremental result must equal a complete $(N+1)$-Gaussian rebuild.

Old-old pair data must be inherited rather than recomputed.

## 5. Incremental pruning

Block-sliced S/H matrices and a subset pair cache must equal a full rebuild of the
surviving basis.

## 6. Cost-aware ranking

The cost-aware layer is allowed to rerank only candidates already admitted by the TDSE
residual shortlist.

It must never promote a candidate below the minimum residual-capture requirement.

Synthetic tests verify that a slightly lower-capture but far better-conditioned
candidate can outrank an ill-conditioned candidate when the conditioning penalty is
activated.

## 7. Release cost gate

The release requires

```text
minimum cost-aware utility = 0.15
```

The accepted event has

```text
utility = 0.24210065440802597
```

and therefore passes the gate.

## 8. Physics regression against v0.14

Because the release benchmark selects the same physical candidate as v0.14, the
performance optimization must preserve the benchmark observables.

Maximum difference across the stored acceptance metrics:

```text
8.412825991399586e-12
```

Release threshold:

```text
1e-9
```

## 9. Factorization reduction

The release requires at least 84% reduction in propagation pair-factorization
equivalents relative to the v0.14 helper structure.

Measured:

```text
84.797 %
```

## 10. Cache reuse

Required cache reuse fraction:

```text
>= 0.60
```

Measured:

```text
0.6408125577100646
```

## 11. Incremental expansion

The accepted candidate's exact conditioning cache must be reused by the matrix
expansion.

Required new factorizations during expansion:

```text
0
```

Measured:

```text
0
```

## 12. Physical accuracy thresholds

```text
initial reduced-density error       <= 0.035
projected-state dynamics error      <= 0.003
target full-density error           <= 0.035
target population L2 error          <= 0.03
coherence phase error               <= 0.0035 rad
generalized norm drift              <= 1e-4
maximum overlap condition number    <= 5e3
```

These remain regression criteria for the analytic LVC benchmark.

## 13. Timing policy

Wall-clock speedup is reported but is not an acceptance criterion.

The portable complexity claims are:

```text
pair-cache algebraic equality
exact factorization-equivalent counts
cache hit/reverse-view counts
incremental add/remove invariants
Big-O operation models
```

## 14. PySCF scope

v0.15 does not claim that the pair cache makes a full molecular PySCF/AIMS calculation
cheap.

The analytic Gaussian layer is optimized here. Electronic-structure costs remain a
separate future budget.

See `V15_PYSCF_COST_BRIDGE.md`.
