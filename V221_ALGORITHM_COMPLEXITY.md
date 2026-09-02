# v0.22.1 algorithmic complexity

Let `s` be the state count, `d` the nuclear coordinate count, `n_h` the number of
finite-difference steps, `N_x` the one-dimensional grid size, `N_t` the step count,
`N` the Gaussian count, and `E` the active Gaussian-pair edge count.

## Component-resolved provider audit

The audit evaluates both spin-free and SOC operators at `2 d n_h` neighboring
geometries. Frame transport and dense residuals cost at most

\[
O(n_h d s^3)
\]

time and `O(d s^2)` working storage when provider cost is excluded. The default
`n_h=3` is fixed. Unlike v0.22.0, this cost covers every coordinate and both components,
not a hard-coded one-dimensional/four-state path.

## Symmetry admission

Unitarity, time-reversal square, and projector algebra use dense products and cost
`O(s^3)` time. Storing `J` and the complete projector family requires `O(s^2)` to
`O(n_P s^2)` space, where `n_P` is the number of physical projectors.

## Exact-grid oracle

Static provider evaluation stores `N_x` dense `s`-state potentials. Potential
eigendecompositions and half-step exponentials are computed once:

\[
O(N_xs^3)\ \text{setup time},\qquad O(N_xs^2)\ \text{storage}.
\]

Each Strang step then costs

\[
O(N_xs^2+sN_x\log N_x),
\]

so a trajectory costs

\[
O\!\left(N_xs^3+N_t(N_xs^2+sN_x\log N_x)\right).
\]

The v0.22.1 implementation now matches this precomputation model exactly. The direct-
product grid remains a validation oracle, not a production multidimensional solver.

## Gaussian-basis and sparse convergence ladders

The release ladder uses prescribed dense bases of size 1, 3, and 5. General dense block
storage scales as `O(N^2 s^2)`. Sparse block storage scales as `O((N+E)s^2)`. The
threshold ladder repeats the sparse propagation a fixed four times and compares each
trajectory with one dense reference; it changes validation cost but not production
asymptotics.

