# v0.26.0 Release Notes

Release date: 2026-08-24

## Theme

v0.26.0 is the reference-first multidimensional release.  It promotes the validated
v0.25.3 one-dimensional adaptive Gaussian TDVP into two nuclear dimensions and binds
it to an independent exact-grid CI+SOC calculation.

## Added

- Vector packet centres and momenta, `q[I,mu]` and `p[I,mu]`.
- Positive coordinate-diagonal widths `alpha[I,mu]=exp(eta[I,mu])` and real
  coordinate-diagonal chirps `beta[I,mu]`.
- Full positive-definite constant nuclear mass matrices.
- Exact multidimensional complex Gaussian moments through degree four.
- Fully implicit-midpoint multidimensional McLachlan TDVP propagation.
- An independent two-dimensional periodic FFT/Strang matrix-wavefunction solver.
- A two-state CI with a constant complex SOC gap.
- A complete four-state Kramers-doublet CI+SOC model.
- A five-state model containing two CI-forming singlets and all three triplet
  projections.
- Multidimensional residual candidates along every signed coordinate and momentum
  axis.
- Projection-controlled multidimensional spawning, pruning, and merge-to-survivor.
- Population plus metric-condition plus velocity-amplification gates for newborn
  shape activation.
- Independent dense-grid quadrature, one-dimensional reduction, exact-grid
  convergence, time-reversal, packet-permutation, electronic-gauge, and signed-axis
  covariance validation.

## Corrective decisions

Population alone is not sufficient for multidimensional newborn activation.  At a
coefficient-row population near `1e-6`, the extra eight shape directions of a 2D
packet can technically pass the absolute SVD cutoff while producing a retained
condition number near `1e9` and a several-hundred-fold velocity amplification.
v0.26.0 therefore requires all three conditions:

1. coefficient-row population at least `1e-6`;
2. retained active metric condition number no larger than `1e8`;
3. active velocity norm no larger than 100 times the dormant-system value.

An additional spawn is not permitted while any earlier newborn remains in the
coefficient-only stage.

## Validation result

- Inherited v0.25.3 gates: 715
- New scientific-validation gates: 80
- New adversarial/core gates: 30
- New gates: 110
- Cumulative gates: **825/825**

The adaptive two-packet trajectory reduces the short-time exact-grid wavefunction
error from approximately `1.20e-3` to `6.73e-4` on the frozen CI+SOC benchmark.
This is a controlled benchmark result, not a general molecular-accuracy claim.

## Explicit boundaries

The following remain false:

- full correlated complex width matrices;
- arbitrary coordinate-rotation covariance of anisotropic packets;
- coordinate-dependent electronic frames in the v0.26.0 TDVP solver;
- unrestricted spawning directions or multiple simultaneous events;
- full AIMS branching;
- absorbing-boundary exact-grid propagation;
- live PySCF molecular-SOC trajectories;
- general ab initio SOC-dynamics accuracy.

The recommended next method layer is a full complex symmetric width matrix, after
which arbitrary orthogonal coordinate covariance can be tested honestly.
