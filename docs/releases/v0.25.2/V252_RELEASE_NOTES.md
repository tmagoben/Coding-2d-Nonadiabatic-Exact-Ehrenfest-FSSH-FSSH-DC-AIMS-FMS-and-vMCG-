# v0.25.2 release notes

Release date: 2026-08-24

v0.25.2 is the adaptive-width/chirp variational release. It extends the exact
one-dimensional v0.25.1 metric layer without opening basis spawning, pruning,
multidimensional widths, moving electronic frames, or molecular-SOC trajectories.

## Added

- `ThawedGaussianSpinorStateV252` with a positive width `alpha_I` and real
  quadratic chirp `beta_I` on every packet.
- Logarithmic width coordinates `eta_I=log(alpha_I)`, so nonlinear iterations cannot
  cross through zero width.
- Exact complex unequal-width/chirp cross moments through degree four.
- Adaptive tangent vectors for coefficient real/imaginary parts, center, momentum,
  log-width, and chirp: `P=2*N_g*N_s+4*N_g` real parameters.
- Width-domain gates for minimum/maximum width, maximum chirp, and maximum per-step
  logarithmic-width change.
- Fully bound implicit-midpoint receipts that recompute the adaptive midpoint state,
  metric, RHS, full SVD, velocity, nonlinear residual, norm, energy, and width domain.
- A closed-form harmonic breathing oracle and exact continuous thawed-Gaussian
  reduction for `q`, `p`, `eta`, and `beta`.
- A coherent-state reduction proving that `alpha=m*omega`, `beta=0` reproduces the
  v0.25.1 frozen-width trajectory.
- 70 validation and 25 adversarial/core gates, yielding 95 new and 630 cumulative
  release checks.

## Numerical correction made during validation

The initial HYBR step tolerance was `1e-11`. Under a constant complex electronic
gauge, the solver reached a `2.8e-16` residual but reported a no-progress status
because its requested step tolerance was below the receipt's structural resolution.
v0.25.2 freezes `xtol=1e-10` while retaining the independent `2e-10` nonlinear
residual gate. The transformed solve then reports normal success without weakening
the achieved residual. Both solver success and the recomputed residual remain
mandatory.

## Claim boundary

Validated: one nuclear coordinate, adaptive scalar widths plus chirps, fixed complete
electronic spinors, constant mass, Hermitian quadratic matrix potentials, compatible
SVD null spaces, and implicit midpoint.

Not validated: spawning/pruning, multidimensional or correlated width matrices,
coordinate-dependent electronic frames in this solver, real molecular-SOC
trajectories, and general ab-initio dynamics accuracy.

