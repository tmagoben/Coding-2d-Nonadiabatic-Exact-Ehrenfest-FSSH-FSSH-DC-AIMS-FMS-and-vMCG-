# v0.16 Algorithmic Complexity and Sparse Scaling Audit

## Symbols

| Symbol | Meaning |
|---|---|
| $N$ | Gaussian/TBF count |
| $E$ | active off-diagonal locality edges |
| $M$ | KD-tree spatial candidate pairs |
| $s$ | electronic states per Gaussian |
| $d$ | nuclear dimension |
| $G$ | TDSE-defect grid points |
| $K$ | residual candidate count |

## Core complexity

**Width preprocessing**

Each Gaussian width needs its smallest eigenvalue:

$$
O(Nd^3).
$$

**Safe KD-tree locality query**

Typical bounded-locality behavior:

$$
O(N\log N+M).
$$

Worst case:

$$
O(N^2)
$$

when the global radius spans almost the entire basis.

**Pair-specific upper-bound checks**

$$
O(Md).
$$

**Exact active pair algebra**

With $E$ active off-diagonal edges:

$$
\boxed{
O((N+E)d^3).
}
$$

**Sparse matrix storage**

Nominal block scaling:

$$
\boxed{
O\left[s^2(N+2E)\right].
}
$$

**TDSE defect**

After the local graph and sparse matrices are available:

$$
O(NGs+sG\log G)
$$

for state reconstruction and FFT-grid Hamiltonian work, plus sparse projected solves.

**Candidate residual ranking**

The retained v0.15 vectorized shortlist costs

$$
O\left(
KG(N+s)+N^2K+N^3
\right).
$$

Local-degree and electronic-cost reranking is only

$$
O(K_{\mathrm{short}}Nd)
$$

for a small shortlist.

## Compact CI benchmark

The release basis is highly overlapping.

```text
average off-diagonal sparsity:
5.38 %

propagation pair factorizations:
14973

v0.15 propagation pair factorizations:
15675

pair reduction versus v0.15:
4.48 %
```

This is intentionally modest.

The sparse graph is solving the same highly local packet representation, not a
synthetically separated basis.

## Endpoint sparse matrix audit

```text
relative S error:
0.005191742661052565

relative H error:
0.003962632349871911

omitted off-diagonal edges:
3

maximum omitted overlap:
0.02057350476995086

maximum omitted H-block Frobenius norm:
0.2318916964441307
```

The H audit is necessary because the overlap cutoff alone is not a universal
Hamiltonian error bound.

## Bounded-locality scaling benchmark

| N | Active edges | Edge fraction | Pair reduction | Dense assembly (s) | Sparse assembly (s) | Speedup |
|---:|---:|---:|---:|---:|---:|---:|
| 20 | 37 | 0.19474 | 72.86% | 0.028987 | 0.004917 | 5.90x |
| 40 | 77 | 0.09872 | 85.73% | 0.112917 | 0.008293 | 13.62x |
| 80 | 157 | 0.04968 | 92.69% | 0.523068 | 0.016066 | 32.56x |

The fitted exponents over $N=20,40,80$ are:

```text
active edges:
1.0425836916313385

KD-tree spatial candidates:
1.0425836916313385

actual pair factorizations:
1.027926617366759

dense canonical pair count:
1.9737662900529327
```

This demonstrates approximately linear local graph growth for this **specific**
bounded-locality chain, while dense pair count remains approximately quadratic.

It is not a universal linear-scaling theorem.

## KD-tree screening effect at N=80

```text
all possible off-diagonal pairs:
3160

KD-tree spatial candidates:
157

globally screened without pair-specific bound:
3003

exact pair checks:
157

canonical pair factorizations including diagonals:
237
```

This is the main difference between a genuinely local graph and simply storing a
sparse matrix after an $O(N^2)$ dense pair scan.

## Sparse direct solver caveat

`scipy.sparse.linalg.spsolve` uses sparse direct factorization.

Its time and memory depend on fill-in.

Therefore v0.16 reports:

```text
matrix nnz
graph edge counts
pair-factorization counts
wall times
```

instead of asserting a false universal complexity exponent for all graph topologies.

## Current measured release timing

```text
graph updates:
1.598122 s

sparse S/H assembly:
0.447088 s

sparse T assembly:
0.881072 s

sparse Cayley solves:
0.045959 s

TDSE-defect work:
0.259131 s

candidate residual ranking:
0.183147 s

local/provider cost reranking:
0.006858 s

total adaptive run:
4.391310 s
```

At $N\approx11$, the sparse machinery is not expected to substantially outperform the
optimized dense v0.15 implementation. Sparse overhead is being introduced for the
larger-locality regime.

## Memory roadmap

The pair cache stores dense $d\times d$ cross covariance data only on evaluated local
pairs:

$$
O((N+E)d^2).
$$

The sparse projected matrices scale as

$$
O[s^2(N+2E)].
$$

For high-dimensional molecular applications, structured widths or factored covariance
representations remain an important future optimization.

## Next scaling bottlenecks

Once pair locality is effective, the likely next bottlenecks are:

1. TDSE-defect candidate grids $O(KG)$;
2. sparse direct-solver fill-in;
3. electronic-structure cache misses;
4. high-dimensional dense Gaussian widths.

Those are better v0.17 targets than simply tightening the overlap cutoff.
