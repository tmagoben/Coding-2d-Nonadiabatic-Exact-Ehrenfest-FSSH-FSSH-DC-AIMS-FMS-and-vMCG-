# v0.14 Algorithmic Complexity and Runtime Audit

This document is intentionally explicit because v0.14 combines several individually
nontrivial algorithms.

## Symbols

| Symbol | Meaning |
|---|---|
| $N$ | nuclear Gaussian/TBF count |
| $s$ | electronic states per Gaussian |
| $d$ | nuclear coordinate dimension |
| $G$ | diagnostic grid points |
| $K$ | candidate Gaussians in one defect search |
| $T$ | propagation time steps |
| $m$ | defect-check interval |

For the primary v0.14 benchmark:

```text
initial N = 10
peak N = 11
average N = 10.925
s = 2
d = 2
defect grid = 40 x 40 = 1600 points
peak K = 560
time steps = 120
defect interval = 10
```

---

## Complexity summary

| Algorithm | Leading time complexity | Leading memory |
|---|---|---|
| exact unequal-width S/H pair build | $O(N^2d^3+N^2s^2)$ | $O(s^2N^2)$ |
| moving-basis $T$ matrix | $O(N^2d^3)$ | $O(s^2N^2)$ |
| dense Cayley solve | $O((sN)^3)$ | $O((sN)^2)$ |
| TDSE-defect check | $O(N^2d^3+(sN)^3+NGs+sG\log G)$ | $O(NGs+sG)$ |
| dynamic candidate ranking | $O(KG(N+s)+N^2K+N^3)$ | $O(KG+NG)$ |
| exact leave-one-out pruning audit | $O(N^3+Ns)$ | $O(N^2)$ |
| full propagation | step cost + periodic defect/search overhead | dynamic |

This is a **dense polynomial-time reference implementation**.

It is not designed for hundreds or thousands of simultaneously active TBFs without
additional locality, sparsity, caching, or iterative linear algebra.

---

## Mathematical/software complexity level

The individual layers have different kinds of complexity.

| Layer | Complexity level | Why |
|---|---|---|
| unequal-width Gaussian integrals | high mathematical / moderate software | complex centroids, dense SPD solves, exact moments |
| nonorthogonal moving-basis TDSE | high | generalized metric, moving connection, norm consistency |
| Cayley propagation | moderate-high | implicit dense solve in a time-dependent metric |
| TDSE-defect reconstruction | high | couples analytic Gaussian basis to an independent FFT-grid Hamiltonian |
| residual candidate ranking | high | nonorthogonal projection + large grid contractions + conditioning |
| adaptive controller | high software complexity | hysteresis, cooldowns, basis budget, growth, pruning, emergency conditioning |
| molecular PySCF extension | very high | adds expensive electronic structure and gauge/state-tracking cost |

So v0.14 is no longer a simple pedagogical propagator.

It is still deliberately readable, but it is an advanced adaptive numerical
algorithm.

---

## Exact Hermitian pair-count reduction

The old ordered-pair builder evaluates:

$$
N^2.
$$

The v0.14 half-builder evaluates:

$$
\frac{N(N+1)}{2}.
$$

At $N=11$:

```text
ordered pairs = 121
Hermitian half-build pairs = 66
reduction = 45.45 %
```

Across the actual variable-basis release trajectory:

```text
Hermitian pair evaluations:
7931

ordered-pair equivalent:
14531

actual reduction:
45.420 %
```

This optimization leaves the matrix elements unchanged to regression precision.

---

## Runtime ledger for the release campaign

Total instrumented adaptive-run wall time:

```text
11.289004 s
```

| Component | Seconds | Fraction of total |
|---|---:|---:|
| Exact S/H pair matrices | 4.885241 | 43.27% |
| Moving-basis T matrices | 4.421099 | 39.16% |
| TDSE-defect evaluations | 0.680050 | 6.02% |
| Candidate ranking | 0.161364 | 1.43% |
| Cayley solves | 0.007415 | 0.07% |
| Pruning audits | 0.002369 | 0.02% |
| Other Python/control/trajectory work | 1.131465 | 10.02% |

Two observations matter.

First, at the current tiny electronic dimension

$$
2N\le22,
$$

the dense Cayley solve is not the practical bottleneck despite its cubic asymptotic
scaling.

Second, the two pairwise Gaussian matrix constructions dominate:

```text
S/H exact pair matrices
+
moving-basis T matrices
```

Together they account for roughly

```text
82.44 %
```

of measured wall time.

That immediately identifies a concrete optimization target for a future release:
**cache/reuse pair algebra across S/H/T construction rather than recomputing the same
Gaussian cross quantities independently.**

---

## Why $T$ remains expensive

The S/H half-build can use Hermiticity.

The moving-basis connection cannot, because

$$
T_{ij}
=
\langle g_i|\dot g_j\rangle
$$

is not Hermitian.

Only

$$
\dot S=T+T^\dagger
$$

is constrained.

Therefore v0.14 correctly keeps the full ordered-pair $T$ construction.

This is why the timing audit shows that optimizing S/H alone does not remove the full
pairwise bottleneck.

---

## Candidate-ranking complexity in the actual event

At the enrichment event:

```text
K = 560
G = 1600
N ≈ 10
s = 2
```

The leading contraction proxy

$$
KG(N+s)+N^2K+N^3
$$

is approximately

```text
10,809,000
```

scalar-scale contraction units.

The prepared candidate grid contains:

```text
896,000 complex128 values
≈ 13.67 MiB
```

By comparison, one $22\times22$ complex matrix is only about

```text
0.0074 MiB
```

in raw `complex128` storage.

Thus for the present benchmark, **candidate-grid memory**, not dense quantum-matrix
memory, is the larger temporary array.

---

## Complexity overhead actually paid by adaptation

The release performs:

```text
matrix builds: 122
moving-basis T builds: 120
Cayley solves: 120
defect evaluations: 13
candidate searches: 1
candidate Gaussians scored: 560
enrichments: 1
pruning audits: 11
actual propagated pruning events: 0
```

Only one enrichment was required.

The candidate-ranking wall time was

```text
0.161364 s
```

or about

```text
1.43 %
```

of total runtime.

The defect checks cost more than the candidate ranking itself because each check
reconstructs the state, builds a moving-basis connection, solves the projected
equation, and applies an FFT-grid Hamiltonian.

---

## Effect of doubling basis size

Ignoring changes in $K$ and $G$:

- pairwise matrix work scales approximately as $N^2$, so $N\to2N$ gives about
  $4\times$ work;
- dense solve work scales as $N^3$, so $N\to2N$ gives about $8\times$ work;
- basis-grid reconstruction scales as $N$, so $N\to2N$ gives about $2\times$ work;
- candidate ranking can grow even faster in practice if $K$ itself is proportional to
  $N$.

If the local dynamic dictionary has

$$
K\propto N,
$$

then the dominant grid contraction

$$
KGN
$$

becomes approximately

$$
O(N^2G).
$$

---

## Effect of increasing nuclear dimension

The exact unequal-width Gaussian pair algebra uses dense $d\times d$ solves/inverses.

Therefore its generic cost contains

$$
O(d^3).
$$

The current benchmark has only

$$
d=2.
$$

A direct full-dimensional molecular implementation with dozens of nuclear coordinates
would therefore require substantial structural optimization even before considering
the electronic-structure cost.

Possible future reductions include:

- diagonal or block-diagonal Gaussian width matrices;
- Cholesky factors cached per pair;
- low-rank coordinate subspaces;
- mode-local Gaussian products;
- sparse/local overlap graphs.

---

## PySCF complexity is not included in these timings

The release benchmark uses an analytic LVC provider.

A real PySCF SA-CASSCF direct-dynamics run adds electronic-structure costs that can
dominate all Gaussian linear algebra.

Those costs depend on:

- basis-function count;
- active orbital/electron count;
- number of states;
- SCF/CASSCF iteration counts;
- gradient/NAC evaluation;
- cross-geometry many-electron overlaps.

Therefore the v0.14 ledger should be interpreted as the complexity of the **Gaussian
adaptive dynamics layer**, not the total cost of an ab initio calculation.

---

## Practical scaling conclusion

At the present benchmark size, v0.14 is comfortably small.

At larger Gaussian counts, three bottlenecks become important in sequence:

```text
1. repeated pairwise Gaussian algebra
2. large K x G residual candidate contractions
3. dense O((sN)^3) coefficient solves
```

A scalable future implementation should therefore prioritize:

```text
pair-data caching
local/sparse overlap graphs
candidate pre-screening
iterative or block linear solvers
batched electronic-structure caching
```

before simply increasing `max_basis`.
