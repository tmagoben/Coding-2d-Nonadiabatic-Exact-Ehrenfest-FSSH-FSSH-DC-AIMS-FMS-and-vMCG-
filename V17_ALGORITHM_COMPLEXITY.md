# v0.17 Algorithmic Complexity and Error-Control Audit

v0.17 makes a distinction that is easy to miss:

> sparse storage, sparse active edges, and sparse **construction cost** are three
> different things.

The release measures all three.

## Symbols

| Symbol | Meaning |
|---|---|
| $N$ | Gaussian/TBF count |
| $E$ | active off-diagonal graph edges |
| $M$ | exact S/H/T-scored local pair candidates |
| $s$ | electronic states per Gaussian |
| $d$ | nuclear dimension |
| $G$ | TDSE-defect grid points |
| $K$ | residual candidate count |
| $m$ | dense-audit interval in propagation steps |

## 1. Geometric candidate generation

Width preprocessing requires smallest eigenvalues of $A_i$:

$$
O(Nd^3).
$$

The safe KD-tree query is typically

$$
O(N\log N+M_{\rm spatial})
$$

for bounded local density, followed by pair-specific bound checks

$$
O(M_{\rm spatial}d).
$$

Worst-case geometry can still yield

$$
M_{\rm spatial}=O(N^2).
$$

## 2. Exact S/H/T edge scoring

Each exact local score reuses the Gaussian pair object but still requires the
cross-Gaussian linear algebra for a new pair.

For $M$ exact-scored pairs,

$$
\boxed{
O(Md^3)
}
$$

is the leading pair-algebra cost.

The score then adds only small $2\times2$ block norms and scalar moving-basis terms.

## 3. Sparse matrix assembly

With $E$ active off-diagonal edges, exact sparse S/H assembly is

$$
O((N+E)d^3+s^2(N+2E)).
$$

After pair data are cached, the block insertion itself is approximately

$$
O(s^2(N+2E)).
$$

## 4. Local score-budget promotion

Suppose $D$ edges are tentatively omitted after exact scoring.

Sorting them by score costs

$$
O(D\log D).
$$

The subsequent promotion loop is linear.

This is lower order than pair algebra for the intended regime.

## 5. Sparse Cayley solve

v0.17 continues to use sparse direct factorization.

Its cost depends on graph topology and fill-in, so the release does **not** assert a
universal $O(N^p)$ sparse-solve exponent.

The portable diagnostics are:

```text
matrix nnz
active edges
solve call count
wall time
```

## 6. Periodic dense audit

One full dense S/H audit costs

$$
\boxed{
O(N^2d^3+s^2N^2).
}
$$

Every $m$ steps, the amortized cost is

$$
\boxed{
O\left(
\frac{N^2d^3+s^2N^2}{m}
\right).
}
$$

For v0.17,

```text
audit interval = 20 steps
dense audits = 8
audit pair factorizations = 506
audit time = 0.148051 s
```

These audits are deliberately accepted as temporary correctness overhead.

They should eventually be replaced by cheaper sampled/operator-residual audit
strategies after the approximation hierarchy is better established.

## 7. Compact CI benchmark timing

```text
graph / S-H-T score work:
8.550292 s

sparse S/H assembly:
0.495324 s

sparse T assembly:
1.460446 s

sparse Cayley solves:
0.085747 s

TDSE-defect work:
0.481460 s

candidate ranking:
0.363455 s

dense audit work:
0.148051 s

total adaptive run:
13.552188 s
```

The graph scoring is now more expensive than v0.16 overlap-only locality on this tiny,
highly overlapping basis.

That is expected: v0.17 pays extra work to measure $H$ and $T$ importance and to
validate the sparse approximation online.

The purpose of v0.17 is **controlled sparsity**, not a benchmark-time speed record at
$N\approx11$.

## 8. Bounded-locality scaling

| N | Active edges | Exact score checks | Pair factorizations | Reduction vs dense | Sparse assembly (s) | Dense assembly (s) | Speedup |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 20 | 37 | 70 | 90 | 57.14% | 0.003547 | 0.057389 | 16.18x |
| 40 | 77 | 150 | 190 | 76.83% | 0.006368 | 0.243558 | 38.25x |
| 80 | 157 | 310 | 390 | 87.96% | 0.012810 | 0.949500 | 74.12x |
| 160 | 317 | 630 | 790 | 93.87% | 0.029290 | 3.823627 | 130.54x |

At $N=160$:

```text
active edges:
317

all possible off-diagonal edges:
12720

exact S/H/T score checks:
630

pair factorizations including diagonals:
790

dense canonical pair count:
12880

pair reduction:
93.87 %

dense/sparse assembly speedup:
130.54 x
```

## 9. Fitted construction exponents

Over $N=20,40,80,160$:

```text
active-edge exponent:
1.032449

KD-tree spatial-candidate exponent:
1.055708

exact S/H/T score-check exponent:
1.055708

pair-factorization exponent:
1.043904

dense canonical pair exponent:
1.979810
```

For this bounded-locality chain, the actual decision cost is close to linear, not just
the final matrix storage.

These fits are benchmark-specific.

## 10. Residual-candidate cost remains

The adaptive TBF candidate ranking inherited from v0.15/v0.16 still costs approximately

$$
O(KG(N+s)+N^2K+N^3).
$$

This remains a major target for the next release.

## 11. Current complexity conclusion

v0.17 establishes the following hierarchy:

```text
dense basis:
    O(N^2) pair construction

bounded-local sparse basis:
    approximately O(N) exact pair scoring in the chain benchmark

online correctness bridge:
    periodic O(N^2) dense audits
```

So v0.17 has **not yet achieved an end-to-end asymptotically sparse algorithm** because
the periodic audit is intentionally dense.

That limitation is explicit rather than hidden.
