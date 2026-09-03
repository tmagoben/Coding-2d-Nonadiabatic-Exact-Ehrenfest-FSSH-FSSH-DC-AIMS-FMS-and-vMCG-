# v0.25.3 release notes

Release date: 2026-08-24

v0.25.3 is the controlled adaptive-basis release. It wraps, but does not alter, the
one-dimensional adaptive-width/chirp TDVP kernel certified in v0.25.2.

## Added

- Analytic residual coupling for displaced position/momentum spawn candidates.
- Orthogonalized candidate novelty, enlarged-overlap rank, and condition gates.
- Full-SVD projection receipts for every topology change.
- Stable packet IDs, packet ages, monotone packet serials, and at most one topology
  event per checkpoint.
- Coefficient-only activation for exactly projected zero-amplitude newborn packets.
- Projection-guarded low-population pruning and high-overlap merge-to-survivor.
- Constant electronic-gauge and packet-permutation covariance validation.
- Exact no-event reduction to the v0.25.2 implicit TDVP step.
- 60 validation and 25 adversarial/core gates: 85 new and 715 cumulative checks.

## Correction found during implementation

An exactly projected newborn must have zero coefficient in a linearly independent
enlarged basis. Its `q`, `p`, log-width, and chirp tangents therefore vanish, making
those four coordinates physically undefined. Attempting the unrestricted v0.25.2
metric on the next midpoint can produce an incompatible numerical null-space RHS.
v0.25.3 does not diagonal-load that metric. Instead, all electronic coefficient
coordinates evolve while dormant shape coordinates remain exactly fixed until the
coefficient-row population reaches `1e-6`. The reduced metric and activation
residual are fully receipted.

## Claim boundary

Validated: controlled one-dimensional residual spawning, coefficient-only newborn
activation, SVD projection, pruning, and merge-to-survivor for fixed-frame Hermitian
quadratic complete-spin Hamiltonians.

Not validated: general or multidimensional AIMS branching, full correlated width
matrices, coordinate-dependent electronic frames in this solver, real molecular-SOC
trajectories, or general ab-initio dynamics accuracy.
