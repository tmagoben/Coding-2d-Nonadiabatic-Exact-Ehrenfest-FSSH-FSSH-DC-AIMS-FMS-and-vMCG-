# v0.26.0 Algorithm Complexity

Let

- `D` be nuclear dimension;
- `S` be electronic-state count;
- `G` be Gaussian packet count;
- `P=2GS+4GD` be the real TDVP parameter count;
- `N=Nx*Ny` be the two-dimensional grid size.

## Exact-grid branch

Pointwise diagonalization is precomputed for a time-independent potential:

$$
O(NS^3)\text{ time},\qquad O(NS^2)\text{ memory}.
$$

Each step uses two pointwise matrix applications and `S` two-dimensional FFTs:

$$
O(NS^2+SN\log N)\text{ time},\qquad O(SN+NS^2)\text{ memory}.
$$

Storing every wavefunction snapshot costs `O(n_store*S*N)`; production calls should
increase `store_every` when full states are unnecessary.

## Gaussian matrix branch

Nuclear overlap and Hamiltonian blocks require all packet pairs:

$$
O(G^2(S^2+D^2))
$$

for quadratic fixed-frame models, with small degree-four polynomial dictionaries.

## TDVP metric branch

The dense real metric contains `P^2` elements.  Direct analytic assembly scales
approximately as

$$
O(P^2D+PGS^2D^4),
$$

where the second term represents Hamiltonian projections with small fixed-degree
polynomials.  Full dense SVD costs

$$
O(P^3)\text{ time},\qquad O(P^2)\text{ memory}.
$$

The nonlinear midpoint cost multiplies metric construction by the number of residual
evaluations.  This release prioritizes auditability over large-`G` optimization.

## Candidate branch

There are `4GD` axis candidates.  Naive candidate evaluation recomputes residual
couplings at the current basis and scales approximately as

$$
O(GD\,[G^2S^2+G^3]).
$$

The current implementation is appropriate for small benchmark bases.  A future
release should cache basis residuals, overlap factorizations, and operator
polynomials across candidates.

## Projection branch

Projecting onto `G'` packets costs

$$
O(G'^3+G'GS+G'^2S)
$$

and stores `O(G'^2+G'G)` complex overlap data.

## Practical boundary

v0.26.0 is a rigor-first low-dimensional research implementation.  Sparse tangent
metrics, iterative solvers, batched candidate caching, and distributed exact grids
are not claimed.
