# Frozen-width multi-Gaussian TDVP (v0.25.1)

v0.25.1 introduces the first genuine coupled variational metric in the release
line. The state is a sum of one-dimensional frozen-width Gaussian packets, each
carrying a complete electronic spinor. Real coefficient parts, imaginary coefficient
parts, packet centers, and packet momenta are solved together by McLachlan variation.

The implementation uses exact Gaussian moments, a full-SVD pseudoinverse with
compatible-null-space auditing, and a fully implicit midpoint residual. It validates
complete even and odd analytic SOC models, exact signed reversal, packet-permutation
and constant-unitary electronic covariance, duplicate-packet null directions,
one-packet harmonic reduction, zero SOC, and second-order timestep refinement.

For the complete derivation, solver thresholds, architecture diagram, complexity,
and evidence tables, see the release-root documents:

- `V251_MULTIGAUSSIAN_TDVP.md`
- `V251_METRIC_AND_SOLVER.md`
- `V251_PROGRAM_ARCHITECTURE.md`
- `V251_ALGORITHM_COMPLEXITY.md`
- `V251_VALIDATION.md`

The exact scope remains one coordinate, frozen widths, a fixed electronic frame,
and a Hermitian quadratic matrix potential. Adaptive widths, spawning/pruning,
multidimensional propagation, coordinate-dependent gauges, and real molecular SOC
trajectories remain later milestones.

