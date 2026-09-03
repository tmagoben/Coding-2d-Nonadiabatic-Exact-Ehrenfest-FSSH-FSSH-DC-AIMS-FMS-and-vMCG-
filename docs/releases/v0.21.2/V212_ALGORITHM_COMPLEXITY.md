# v0.21.2 Algorithmic Complexity

Let $N$ be the Gaussian count, $s$ the electronic block dimension, $d$ the nuclear
dimension, $E$ the active Gaussian-edge count, and $M$ the geometrically admitted pair
count.

## Unequal-width pair algebra

Each general width pair requires solves/factorizations involving $d\times d$ matrices,
so the nuclear algebra remains

$$
O(d^3)
$$

per admitted pair. The electronic block algebra adds

$$
O(s^3)
$$

for polar/SVD and dense block products.

Thus the pair layer remains approximately

$$
O\left[M(d^3+s^3+C_{ES})\right].
$$

## Self-consistent guidance

For each Gaussian, the force expectation over $d$ physical operator derivatives costs

$$
O(ds^2)
$$

once the electronic point is available. For $N$ Gaussians,

$$
O(Nds^2)
$$

plus provider cost.

With $k$ predictor/corrector iterations, one nuclear step multiplies the endpoint block
matrix/provider work by approximately $k+1$. v0.21.2 defaults to two corrector
iterations and uses this path as a correctness/hardening reference, not as a claim of
optimal production scheduling.

## Adaptive birth

Zero-block insertion costs

$$
O(Ns)
$$

to allocate/copy the coefficient vector and does no projection solve.

## Block pruning

Deleting one TBF leaves a retained metric of dimension $(N-1)s$. A dense projection
solve is formally

$$
O(((N-1)s)^3).
$$

This is currently a robust reference implementation. Future large-basis work can reuse
factorizations or exploit sparse Schur complements.

## Electronic observables

For an active-edge observable matrix, storage is

$$
O[s^2(N+E)].
$$

Centroid operator assembly costs approximately

$$
O[E(s^3+C_O)],
$$

where $C_O$ is the provider/evaluator cost of the physical observable matrix.

## Full-subspace diagnostics

One $s\times s$ Procrustes SVD costs

$$
O(s^3).
$$

Nearest trusted-anchor lookup retains the buffered KD-tree cost inherited from v0.20.

## Key point before SOC

None of these asymptotic layers depends on whether the electronic Hamiltonian contains
SOC. A later SOC backend changes the electronic operator content and usually $s$, but
not the architecture of the Gaussian pair graph, block storage, metric propagation, or
observable algebra.
