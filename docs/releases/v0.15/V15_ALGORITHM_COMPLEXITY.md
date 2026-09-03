# v0.15 Algorithmic Complexity and Cache Audit

This release treats complexity as part of the numerical method.

## 1. Symbols

| Symbol | Meaning |
|---|---|
| $N$ | nuclear Gaussian/TBF count |
| $s$ | electronic states per Gaussian |
| $d$ | nuclear coordinate dimension |
| $G$ | diagnostic grid points |
| $K$ | dynamic candidates |
| $h$ | cost-estimation horizon in time steps |
| $q$ | defect checks inside the cost horizon |

For the release benchmark:

```text
initial N = 10
peak/final N = 11
average N = 10.925
s = 2
d = 2
G = 40 x 40 = 1600
K at enrichment = 560
residual shortlist = 8
```

## 2. Leading-order complexity

| Component | v0.15 time complexity | Main memory |
|---|---|---|
| one pair snapshot | $O(N^2d^3)$ | $O(N^2d^2)$ |
| cached exact S/H assembly | $O(N^2d^3+N^2s^2)$ | $O(s^2N^2+N^2d^2)$ |
| cached moving-basis T after pair moments | $O(N^2s^2)$ | $O(s^2N^2)$ |
| one-TBF fixed-snapshot append | $O(Nd^3+Ns^2)$ | $O(Nd^2+s^2N)$ incremental |
| fixed-snapshot prune after audit | $O((sN)^2)$ slicing/copying | $O((sN)^2)$ |
| Cayley/Galerkin solve | $O((sN)^3)$ | $O((sN)^2)$ |
| TDSE-defect check | $O(NGs+sG\log G+N^2s^2+(sN)^3)$ after endpoint cache | $O(NGs+sG)$ |
| candidate residual ranking | $O(KG(N+s)+N^2K+N^3)$ | $O(KG+NG)$ |
| exact pruning audit | $O(N^3+Ns)$ | $O(N^2)$ |

The implementation remains a dense polynomial-time reference algorithm.

## 3. Why the pair cache changes the constant dramatically

v0.14 used several helper routines that each independently solved or inverted

$$
A_i+A_j.
$$

For one canonical S/H pair the factorization-equivalent bookkeeping count was 7.

v0.15 uses one multi-right-hand-side solve per canonical pair.

For the moving-basis matrix, v0.14 used 3 pair operations for every ordered pair.
v0.15 uses one canonical pair solve and exact conjugate reversal for the opposite
orientation.

The propagation ledger reports:

```text
v0.14 factorization-equivalent baseline:
103103

v0.15 propagation pair factorizations:
15675

avoided:
87428

reduction:
84.80 %
```

The candidate shortlist used another

```text
88
```

factorizations for exact conditioning. They are kept separate from propagation work.

## 4. Pair-cache reuse

```text
pair requests: 23826
direct cache hits: 8063
reverse-orientation views: 7205
inherited pairs reused: 440
cache reuse fraction: 64.08 %
```

The cache reuse fraction counts direct hits plus reverse views divided by all pair
requests. It is not the same as the factorization-reduction fraction.

## 5. Incremental growth

At a fixed endpoint, adding one Gaussian requires only

$$
N+1
$$

new canonical child pairs.

But the release goes further: exact conditioning of the residual shortlist already
creates an expanded cache for each shortlisted candidate. The accepted candidate keeps
that cache.

Therefore the release event required

```text
0
```

new factorizations during the actual matrix expansion.

That is an explicit acceptance criterion.

## 6. Incremental pruning

Pruning keeps the same surviving basis functions at the same geometry, so all
surviving pair data remain valid.

After the $O(N^3)$ leave-one-out pruning audit:

```text
S -> block slice
H -> block slice
Snuc -> block slice
pair cache -> subset/remap
```

No pair integral is recomputed.

## 7. Cost-aware growth complexity

For a one-TBF addition over $h$ future steps,

$$
\Delta F_P=2h(N+1)
$$

additional endpoint+midpoint pair factorizations are expected.

The relative pair overhead is

$$
r_P=2/N.
$$

The dense solve overhead is

$$
r_C=
rac-159.
$$

The release combines pair, solve, defect-check, and conditioning contributions into a
normalized cost $C_c$ and chooses by

$$
U_c=f_c/C_c.
$$

The cost reranking itself is only $O(K_{short})$ after the expensive residual
shortlist has been generated.

In the release:

```text
residual shortlist = 8
cost-aware utility = 0.24210065440802597
normalized incremental cost = 0.6223729662698816
estimated incremental horizon seconds = 0.03330978168687729
```

## 8. Dense-solve scaling

The current peak electronic dimension is only

$$
sN=22.
$$

A generic dense solve is $O((sN)^3)$, but at this size its absolute cost is still
small.

The ledger records:

```text
Cayley calls: 120
Cayley cubic units: 1251280
defect projected solves: 13
defect cubic units: 135776
```

As $N$ increases, this cubic term will eventually overtake the pair algebra.

## 9. Measured timing breakdown

Total adaptive-run wall time in this build:

```text
4.207306 s
```

| Component | Seconds | Fraction |
|---|---:|---:|
| endpoint S/H cached builds | 1.133081 | 26.93% |
| moving-basis T builds | 1.668603 | 39.66% |
| TDSE defect evaluations | 0.254956 | 6.06% |
| candidate residual ranking | 0.156593 | 3.72% |
| cost-aware reranking | 0.000133 | 0.00% |
| Cayley solves | 0.006800 | 0.16% |
| pruning audits | 0.002573 | 0.06% |
| other control/trajectory/grid work | 0.984568 | 23.40% |

Wall time depends on hardware, Python/BLAS libraries, and system load. It is diagnostic,
not an acceptance threshold.

## 10. v0.14 versus v0.15 timing

Saved benchmark timings:

```text
v0.14 adaptive run: 11.289004 s
v0.15 adaptive run: 4.207306 s

diagnostic speedup: 2.683 x
runtime reduction: 62.73 %
```

The key portability result is not the wall-clock ratio. It is that the physical result
is unchanged while repeated pair linear algebra is removed exactly.

## 11. Memory scaling of the pair cache

Each canonical pair stores leading mathematical data of size

$$
O(d^2)
$$

because the covariance matrix is dense.

Thus pair-cache memory is

$$
O(N^2d^2).
$$

At the present $d=2$, this is tiny.

For many nuclear dimensions, dense covariance storage can become significant. Future
high-dimensional versions should consider structured widths, factored covariances, or
local mode blocks.

## 12. Candidate-grid memory

Residual ranking stores $K$ candidate Gaussians on $G$ grid points:

$$
O(KG).
$$

For the release event,

```text
K = 560
G = 1600
K G = 896000 complex values
```

At `complex128`, the raw candidate grid is approximately

```text
13.67 MiB
```

before temporary matrix products.

## 13. Scaling roadmap

The likely bottleneck sequence as the basis grows is:

```text
small N:
    pair algebra / Python pair loops

moderate N:
    candidate K x G contractions + pair cache memory

larger N:
    dense O((sN)^3) coefficient solves

ab initio:
    electronic-structure evaluations can dominate all of the above
```

The next scalable developments should therefore be:

```text
local/sparse overlap graphs
persistent pair-cache invalidation by locality
batched pair algebra
structured width matrices
iterative/block linear solves
electronic-structure cache accounting
```

rather than simply raising `max_basis`.
