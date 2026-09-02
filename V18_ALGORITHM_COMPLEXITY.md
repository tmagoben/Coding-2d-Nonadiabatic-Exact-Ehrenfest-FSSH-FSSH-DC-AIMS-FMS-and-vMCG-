# v0.18 Algorithmic Complexity

v0.18 targets two specific costs left by v0.17:

1. repeated dense audit work;
2. peak $K\times G$ candidate-grid memory.

It does not claim that every remaining algorithmic component is asymptotically sparse.

## Symbols

| Symbol | Meaning |
|---|---|
| $N$ | Gaussian basis size |
| $E$ | active sparse graph edges |
| $M$ | local exact S/H/T-scored pairs |
| $K$ | adaptive Gaussian candidates |
| $B$ | candidate batch size |
| $G$ | diagnostic grid points |
| $d$ | nuclear dimension |
| $J$ | sampled omitted edges |
| $s$ | electronic states |

## Sparse graph construction

The v0.17 structure remains:

$$
O(Nd^3+N\log N+Md)
$$

for width preprocessing, KD-tree search, and geometric bounds, followed by roughly

$$
O(Md^3)
$$

exact pair algebra.

Worst-case dense configurations remain quadratic in $N$.

## Sparse projected matrices

For $E$ active off-diagonal edges,

$$
O((N+E)d^3)
$$

pair algebra is required, with block storage approximately

$$
O[s^2(N+2E)].
$$

Sparse direct-solver fill remains topology dependent.

## Batched candidate ranking

The arithmetic work remains approximately

$$
O(KG(N+s)+N^2K+N^3).
$$

The key v0.18 change is memory.

Earlier peak candidate-grid storage:

$$
O(KG).
$$

v0.18 peak candidate-grid storage:

$$
\boxed{
O(BG).
}
$$

Canonical measured values:

```text
dense candidate-grid elements:
1044800

batched peak elements:
25600

peak reduction:
97.55 %
```

This is a memory-complexity improvement, not a reduction in total grid contractions.

## Sampled audit cost

A normal v0.18 sparse audit evaluates only $J$ omitted edges:

$$
\boxed{
O(Jd^3).
}
$$

The priority search uses the existing geometric locality structure.

In the bounded-locality audit benchmark, exactly 16 pairs are sampled for each basis
size:

| $N$ | All off-diagonal pairs | Sampled pairs | Sample fraction |
|---:|---:|---:|---:|
| 20 | 190 | 16 | 8.4211% |
| 40 | 780 | 16 | 2.0513% |
| 80 | 3160 | 16 | 0.5063% |
| 160 | 12720 | 16 | 0.1258% |

At $N=160$, the audit touches only

$$
0.1258\%
$$

of all possible off-diagonal pairs.

This does not include the cost of constructing the sparse graph itself; it isolates
normal **audit** scaling.

## Dense sentinel cost

A complete sentinel remains

$$
O(N^2d^3+s^2N^2).
$$

v0.18 uses exactly two release sentinels.

Canonical:

```text
v0.17 dense audit pair factorizations:
506

v0.18 dense sentinel pair factorizations:
146

reduction:
71.15 %
```

The end-to-end algorithm is therefore still not mathematically free of $O(N^2)$
validation work, but that work is no longer incurred every audit interval.

## Full-wavefunction metrics

Reconstructing the spinor Gaussian state on $G$ points costs approximately

$$
O(NGs).
$$

Density and moment diagnostics are then $O(Gs)$.

The v0.18 release evaluates these at stored convergence checkpoints rather than every
internal propagation step.

## Canonical timing ledger

```text
wall time for independent canonical coordinate:
9.741 s

internal sparse runner total:
9.740 s

graph/S-H-T score work:
5.302 s

candidate ranking:
0.483 s

sampled audits:
0.0162 s

dense sentinel audits:
0.0224 s
```

Wall times are diagnostic only.

## Remaining v0.19-era scaling bottlenecks

After v0.18 the major computational issues are:

1. total $O(KG)$ residual-candidate arithmetic despite bounded peak memory;
2. sparse direct-solver fill-in;
3. high-dimensional dense Gaussian width algebra;
4. eventual electronic-structure provider cost;
5. sampled-audit statistical/calibration strategy for molecular systems.

These should be addressed before introducing additional complex physics such as SOC.
