# v0.27.0 Release Notes

Release date: 2026-08-25

## Theme

v0.27.0 is the full-correlated-width release. It replaces the coordinate-diagonal
Gaussian shape manifold in v0.26.0 with full symmetric width and chirp matrices,
while preserving the real McLachlan principle, full-SVD metric solution, and fully
implicit midpoint propagation.

## Added

- Full positive-definite width matrices `Gamma[I,:,:]=exp(E[I,:,:])`.
- Full real symmetric quadratic chirp matrices `B[I,:,:]`.
- Frobenius-orthonormal symmetric-matrix coordinates (`svec`).
- Exact Frechet derivatives of the symmetric matrix exponential.
- Exact multivariate complex-normal Wick moments through total degree four.
- Analytic correlated overlap, kinetic, quadratic-potential, tangent-metric, and
  residual-coupling elements.
- Arbitrary orthogonal coordinate covariance, including proper rotations and
  reflections, for states, models, matrices, metrics, velocities, and midpoint steps.
- Correlated spawn, projection, merge, prune, and newborn-activation procedures.
- Intrinsic principal-axis candidates for nondegenerate widths, with deliberate
  fail-closed behavior for degenerate or near-degenerate eigenspaces.
- Independent dense FFT quadrature and exact matrix-Riccati validation oracles.

## Corrective decision

An inactive packet is frozen in log-width coordinates during the nonlinear solve.
Reconstructing `exp(log(Gamma))` can differ from the original matrix by roundoff even
when the coordinate is unchanged. The endpoint now restores dormant `q`, `p`,
`Gamma`, and `B` directly from the source state, making the exact-freeze contract
bitwise true.

## Validation result

- Inherited v0.26.0 gates: 825
- New scientific-validation gates: 100
- New adversarial/core gates: 35
- New gates: 135
- Cumulative gates: **960/960**

The matrix-Riccati oracle generates nonzero off-diagonal shape velocity in a rotated
harmonic Hamiltonian. The full manifold reproduces that velocity to floating-point
precision and its implicit-midpoint trajectory converges at second order. A
coordinate-diagonal manifold cannot represent this correlation rate.

## Explicit boundaries

The following remain false:

- coordinate-dependent electronic-frame covariance in this TDVP kernel;
- live molecular-SOC evaluation inside a trajectory;
- general ab initio SOC-dynamics accuracy;
- optimized directions inside a degenerate width eigenspace;
- simultaneous multiple topology events;
- full AIMS branching semantics.

The recommended next iteration is a coordinate-dependent electronic-frame/connection
layer, initially against analytic moving-gauge oracles, before admitting live
molecular SOC into the variational trajectory.
