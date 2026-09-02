# v0.26.0 Multidimensional Lifecycle Policy

## Checkpoint order

Exactly one of the following may occur:

1. merge;
2. prune;
3. spawn;
4. no event.

The first accepted event ends the checkpoint.

## Spawn candidates

For every packet `I` and coordinate `a`, four candidates are generated:

$$
q_{Ia}\rightarrow q_{Ia}\pm\alpha_{Ia}^{-1/2},
\qquad
p_{Ia}\rightarrow p_{Ia}\pm\alpha_{Ia}^{1/2}.
$$

All untouched coordinates, widths, and chirps are inherited from the source packet.
Candidates are evaluated in a deterministic canonical order.

Let `b_I=<g_I|g_c>` and solve `S a=b` by full SVD.  The novelty is

$$
\nu=1-b^\dagger S^+b.
$$

For `R=dPsi/dt+i H Psi`, the orthogonalized coupling is

$$
r_\perp=\frac{\langle g_c|R\rangle
-a^\dagger\langle\mathbf g|R\rangle}{\sqrt\nu}.
$$

The score is `||r_perp||`.  Admission requires sufficient score and novelty, full
rank of the enlarged overlap, and its condition number below the configured gate.

## Projection

For target basis `h`,

$$
S^{(h)}C'=B C,
\qquad B_{IJ}=\langle h_I|g_J\rangle.
$$

The full-SVD solution records target rank and condition, source and projected norm,
source/projected overlap, normalized fidelity, relative projection loss, and energy
jump.  Reduced-basis events are normalized only after these pre-normalization
quantities are recorded.

An enlarged spawn contains the entire old basis, so the old coefficients are copied
exactly and the newborn row is exactly zero.  Projection loss and energy jump are
therefore exactly zero in the released representation.

## Merge

Pairs are examined deterministically.  A pair must exceed overlap `0.997`.  The
smaller coefficient row is removed, and the full state is projected onto the
survivor basis.  Averaged packet geometry is not used.

## Prune

A packet must be at least 64 steps old and have coefficient-row population below
`1e-10`.  Removal is still rejected unless projection loss and energy jump pass.

## Newborn activation

All coefficient real/imaginary directions are active immediately.  Shape directions
remain zero-velocity until:

- coefficient-row population is at least `1e-6`;
- the trial active metric condition is at most `1e8`;
- trial velocity norm is at most 100 times the dormant-system norm.

Previously activated packet shapes remain active.  Candidate packets are tested one
at a time, and no additional spawn is permitted while any newborn is dormant.

## Closed cases

- more than one event per checkpoint;
- compound split/merge operations;
- arbitrary residual-optimized direction searches;
- energy-adjusted child momenta;
- full AIMS parent/child branching semantics;
- simultaneous activation of several marginal shape blocks.
