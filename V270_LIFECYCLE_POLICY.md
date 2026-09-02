# v0.27.0 Correlated Basis Lifecycle Policy

## Event order

At an adaptation checkpoint, at most one topology change is allowed in the frozen
order:

1. merge;
2. prune;
3. residual-driven spawn;
4. otherwise record no event.

## Spawn candidates

For each packet with nondegenerate width eigenpairs $(\lambda_a,u_a)$, candidates
are displaced along every signed intrinsic axis:

$$
q' = q \pm \frac{c_q}{\sqrt{\lambda_a}}u_a,
\qquad
p' = p \pm c_p\sqrt{\lambda_a}u_a.
$$

The full width and chirp matrices are copied. A candidate must pass novelty,
enlarged-overlap rank, condition number, and orthogonalized TDVP-residual capture
gates. The best admitted residual score wins; stable IDs and ages are deterministic.

## Projection receipts

Merge and prune solve the fixed-time full-SVD least-squares projection. Admission
requires bounded relative projection loss and energy jump. Receipts include norm,
fidelity, source/projected overlap, spectrum, cutoff, rank, condition number, and
energy before/after.

## Newborn activation

An exactly projected newborn has a zero coefficient. Its electronic coefficients
are active immediately, but its complete shape block

$$
(q,p,\operatorname{svec}E,\operatorname{svec}B)
$$

remains dormant until population, retained metric condition, and velocity
amplification gates pass. No additional spawn is allowed while an earlier newborn
is dormant. Frozen endpoint matrices are restored bitwise from the source state.

## Closed procedures

- no laboratory-axis fallback for degenerate widths;
- no multiple simultaneous events;
- no unprojected deletion or merge;
- no diagonalization of a correlated packet to fit an older state type;
- no claim of full AIMS spawning/branching semantics.
