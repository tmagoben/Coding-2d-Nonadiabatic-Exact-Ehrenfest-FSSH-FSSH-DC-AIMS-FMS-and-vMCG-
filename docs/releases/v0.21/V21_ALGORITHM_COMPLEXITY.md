# Algorithmic complexity

Let:

- $N$: nuclear Gaussian count;
- $s$: electronic block dimension;
- $E$: active off-diagonal Gaussian edges;
- $M$: geometrically admitted candidate pairs;
- $d$: nuclear dimension;
- $N_c$: electronic-cache size;
- $C_{\mathrm{ES}}$: cost of one electronic-structure point.

## 1. Electronic operator transformations

For an $s\times s$ operator,

$$
G^\dagger A G
$$

costs

$$
O(s^3).
$$

For $d$ physical derivative operators and $d$ connection matrices, the local
gauge-transform work is

$$
O(ds^3).
$$

## 2. Full-subspace tracking

One SVD of an $s\times s$ overlap matrix costs

$$
\boxed{
O(s^3).
}
$$

This is preferable to factorial root permutation search when the entire manifold is
treated as one subspace.

## 3. One block Gaussian pair

General unequal-width Gaussian algebra is typically $O(d^3)$.

Electronic polar decomposition and dense block multiplications are $O(s^3)$.

Thus one admitted pair costs approximately

$$
\boxed{
O(d^3+s^3+C_{\mathrm{ES}}).
}
$$

For $M$ admitted pairs:

$$
O\left[M(d^3+s^3+C_{\mathrm{ES}})\right].
$$

## 4. Block sparse storage

Each active Gaussian edge carries an $s\times s$ block.

Therefore storage scales as

$$
\boxed{
O[s^2(N+E)].
}
$$

The release state-dimension benchmark confirms this exactly for the Hamiltonian blocks:

| s | Total dimension | Active edges | H nnz | H density |
|---:|---:|---:|---:|---:|
| 2 | 48 | 45 | 456 | 0.197917 |
| 4 | 96 | 45 | 1824 | 0.197917 |
| 8 | 192 | 45 | 7296 | 0.197917 |

The measured values of `H_nnz / s^2` are:

```text
[114.0, 114.0, 114.0]
```

## 5. Sparse pair search

For bounded spatial locality,

$$
M=O(N)
$$

is possible.

Worst-case strongly overlapping Gaussian ensembles remain

$$
M=O(N^2).
$$

The framework does not hide this worst case.

## 6. Coefficient propagation

The coefficient dimension is

$$
M_C=Ns.
$$

A dense direct solve would scale as

$$
O((Ns)^3)
$$

time and

$$
O((Ns)^2)
$$

memory.

The implementation uses sparse matrices; actual sparse factorization complexity depends
on graph topology and fill-in, so no universal sparse exponent is claimed.

## 7. Indexed electronic cache

Exact geometry hits are hash lookups.

Nearest trusted-anchor queries use the buffered KD-tree inherited from v0.20:

$$
O(\log N_c+B d)
$$

between tree rebuilds.

## 8. Dynamic topology cost

The release crossing test uses 51 frames and records:

```text
exact pair checks:
749

entered edges:
15

exited edges:
9
```

This validates repeated graph mutation rather than a single static sparse assembly.

## 9. Dominant remaining costs

For realistic molecular calculations, likely dominant layers are:

1. SCF/SA-CASSCF/gradient/NAC electronic calculations;
2. many-electron cross-geometry overlap;
3. sparse factorization fill-in;
4. adaptive Gaussian candidate generation;
5. complete production-level matrix elements if full AIMS is pursued.

Optional future SOC should alter the electronic operator content, not the asymptotic
structure of the Gaussian engine.
