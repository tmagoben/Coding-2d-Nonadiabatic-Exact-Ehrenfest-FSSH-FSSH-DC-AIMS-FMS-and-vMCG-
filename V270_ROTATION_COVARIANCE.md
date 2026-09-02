# v0.27.0 Orthogonal Coordinate Covariance

## Frozen convention

Coordinates are row vectors related by

$$
R_{\rm old}=R_{\rm new}O,
\qquad OO^T=I.
$$

State and model data transform as

$$
q'=qO^T,\quad p'=pO^T,\quad
E'=OEO^T,\quad B'=OBO^T,
$$

$$
M'=OMO^T,\quad H'_1{}_a=\sum_bO_{ab}H_{1,b},
\quad H'_2{}_{ab}=\sum_{cd}O_{ac}O_{bd}H_{2,cd}.
$$

This convention supports every orthogonal $O$, including $\det O=-1$ reflections.

## Packed tangent map

Let $L(O)$ contain identity blocks for electronic coefficients, $O$ blocks for every
center and momentum, and the induced `svec` representation for every $E$ and $B$.
Because `svec` is Frobenius orthonormal,

$$
L^TL=I.
$$

The validated covariance relations are

$$
A'=LAL^T,\qquad b'=Lb,\qquad \dot\theta'=L\dot\theta.
$$

The implicit midpoint endpoint obeys the same map, not merely its infinitesimal
velocity. Overlap and Hamiltonian matrices, energy, and reduced electronic density
are invariant when the state and model are transformed together.

## Validation coverage

v0.27.0 directly checks:

- `svec` rotation and reflection matrices are orthogonal;
- matrix logarithm, exponential, and Frechet derivative commute with rotations;
- cross overlaps and full covariance moments are invariant;
- wavefunction values agree at corresponding physical coordinates;
- model, overlap, Hamiltonian, energy, and density transform correctly;
- metric, forcing, SVD velocity, and implicit midpoint endpoint are covariant;
- packet permutations and constant complex electronic gauges remain covariant;
- nondegenerate principal-axis candidate scores and the selected spawn event rotate.

## Degenerate eigenvectors

A degenerate width eigenspace has no unique principal axes. Selecting an eigensolver's
arbitrary basis would break deterministic physical covariance. Candidate generation
therefore returns no principal-axis candidates for a packet when adjacent width
eigenvalues differ by less than the frozen relative-gap threshold. Optimizing a
direction inside such an eigenspace is a future method, not a hidden fallback.
