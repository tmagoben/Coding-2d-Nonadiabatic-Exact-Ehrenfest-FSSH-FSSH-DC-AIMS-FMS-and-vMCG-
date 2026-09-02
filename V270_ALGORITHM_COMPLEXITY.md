# v0.27.0 Algorithm Complexity

Let

- $D$ be nuclear dimension;
- $S$ be electronic-state count;
- $G$ be Gaussian packet count;
- $K=D(D+1)/2$ be the number of symmetric-matrix coordinates;
- $P=2GS+2GD+2GK=G(2S+D^2+3D)$ be the real TDVP parameter count;
- $M_4={D+4\choose4}$ be the number of raw moments through degree four.

## Pair algebra

Each packet pair requires a dense complex $D\times D$ solve and moments through a
fixed total degree:

$$
O(D^3+M_4D)
$$

time and $O(D^2+M_4)$ temporary memory. All overlap/Hamiltonian blocks therefore
scale as

$$
O\!\left(G^2[D^3+M_4D+S^2M_4]\right).
$$

## Tangent metric

The dense real metric stores $P^2$ values. With cached pair moments and operator
polynomials, direct assembly is approximately

$$
O(P^2M_4+PGS^2M_4)
$$

before the dense full SVD,

$$
O(P^3)\text{ time},\qquad O(P^2)\text{ memory}.
$$

The nonlinear midpoint multiplies this work by the number of residual evaluations.
v0.27.0 is an audit-first small-basis implementation; it does not claim sparse or
matrix-free scaling.

## Matrix exponential coordinates

Each width log/exp or Frechet derivative costs $O(D^3)$. Forming all $E$ tangents
costs $O(GKD^3)$ without shared divided-difference caching. This is negligible in
the current 2D validation but becomes material at large $D$.

## Candidate lifecycle

Nondegenerate principal-axis generation produces $4GD$ candidates. Naive residual
evaluation scales roughly as

$$
O\!\left(GD\,[G^2S^2M_4+G^3]\right),
$$

because each candidate uses residual projections plus overlap conditioning. A
projection onto $G'$ packets costs $O(G'^3+G'GS+G'^2S)$.

## Optimization direction

The next performance layer should cache pair moment tables, tangent polynomials,
basis residuals, overlap factorizations, and matrix-exponential divided differences.
Sparse/iterative metric solution should be introduced only with new equivalence and
null-space gates; it is not silently substituted in this release.
