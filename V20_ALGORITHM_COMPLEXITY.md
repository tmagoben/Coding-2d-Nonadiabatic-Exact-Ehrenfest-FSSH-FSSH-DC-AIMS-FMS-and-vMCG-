# v0.20 Algorithmic Complexity

## Symbols

| Symbol | Meaning |
|---|---|
| $N$ | Gaussian/TBF count |
| $P=N(N-1)/2$ | dense pair count |
| $M$ | geometrically admitted molecular pair candidates |
| $E$ | active molecular edges |
| $J$ | sampled omitted edges per audit |
| $N_c$ | trusted electronic-cache size |
| $B$ | recent KD-tree buffer size |
| $C_{\mathrm{ES}}$ | cost of one electronic-structure point |

## Candidate discovery

KD-tree discovery is typically

$$
O(N\log N+M_{\mathrm{spatial}})
$$

under bounded local density, followed by pair-specific geometric bounds.

Worst-case highly overlapping ensembles remain $O(N^2)$.

## Molecular pair scoring

Every admitted candidate may require one molecular pair-centroid electronic snapshot:

$$
\boxed{
O(MC_{\mathrm{ES}}).
}
$$

The measured bounded-locality pair-check exponent is

```text
1.0324491201728092
```

rather than the dense pair exponent

```text
2.0213244854837162.
```

## Sparse matrix storage

Only diagonal entries and $E$ off-diagonal edges are stored:

$$
\boxed{
O(N+E)
}
$$

per scalar $S/H/T$ matrix.

The measured active-edge exponent is

```text
1.117285787300988.
```

## Sampled molecular audit

With fixed sample size $J$,

$$
O(JC_{\mathrm{ES}})
$$

new electronic work is required in the worst case per audit.

The canonical release samples at most eight omitted pairs at each checkpoint.

## Dense sentinel

A full validation sentinel remains

$$
O(PC_{\mathrm{ES}}).
$$

It is intentionally isolated in a separate provider/cache and is not the production
sparse path.

## Indexed tracking cache

Nearest-anchor lookup is approximately

$$
O(\log N_c+B n_q)
$$

between buffered KD-tree rebuilds.

Rebuilds cost approximately $O(N_c\log N_c)$ when triggered.

Canonical diagnostics:

```text
nearest queries: 2148
KD queries: 2148
buffer checks: 16086
rebuilds: 135
```

## Bounded-locality scaling

| N | Active edges | Molecular pair checks | Dense pairs | Pair reduction | New electronic points |
|---:|---:|---:|---:|---:|---:|
| 20 | 15 | 37 | 190 | 80.53% | 57 |
| 40 | 32 | 77 | 780 | 90.13% | 117 |
| 80 | 71 | 157 | 3160 | 95.03% | 237 |
| 160 | 152 | 317 | 12720 | 97.51% | 477 |

Fitted exponents:

```text
active edges: 1.117285787300988
pair checks: 1.0324491201728092
new electronic points: 1.0213244854837158
dense pairs: 2.0213244854837162
```

At $N=160$:

```text
pair checks: 317
dense pairs: 12720
pair reduction: 97.51 %
matrix sparsity: 98.81 %
```

## Canonical propagation work

```text
sparse provider cache misses: 2149
dense provider cache misses: 4239
electronic-point reduction: 49.30 %

diagnostic sparse wall time: 3.945 s
diagnostic dense wall time: 7.105 s
diagnostic wall speedup: 1.80 x
```

The electronic-point reduction is the portable metric. Wall time is environment
specific.

## Remaining scaling issues

The former all-pairs molecular centroid construction is no longer the dominant
architectural blocker. Remaining high-cost layers are:

1. actual SCF/SA-CASSCF/gradient/NAC cost;
2. many-electron cross-geometry overlap cost;
3. sparse-factorization fill-in;
4. KD-tree rebuild overhead for extremely large caches;
5. complete AIMS matrix elements if pursued.
