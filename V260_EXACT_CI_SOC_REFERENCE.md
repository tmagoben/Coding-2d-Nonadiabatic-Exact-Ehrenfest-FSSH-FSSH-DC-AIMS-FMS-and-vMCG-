# v0.26.0 Exact CI+SOC Reference

## Purpose

The exact-grid implementation is an independent oracle.  It imports no Gaussian
moment, TDVP metric, implicit solver, projection, or spawning routine.  Agreement
therefore cannot be produced by two paths sharing the same variational algebra.

## Split operator

For a time-independent matrix potential,

$$
e^{-i\hat H\Delta t}=
e^{-i\mathbf V\Delta t/2}
e^{-i\hat T\Delta t}
e^{-i\mathbf V\Delta t/2}+O(\Delta t^3).
$$

At each grid point, `V` is diagonalized with a Hermitian eigensolver.  The kinetic
phase on the FFT wavevector grid is

$$
\exp\left[-\frac{i\Delta t}{2}\mathbf k^TM^{-1}\mathbf k\right].
$$

Both substeps are unitary.  The released boundary condition is periodic; edge-strip
probability is therefore an acceptance gate.

## Two-state CI+SOC model

$$
\mathbf V_{2}(x,y)=
\frac12(\omega_x^2x^2+\omega_y^2y^2)\mathbf 1
+\kappa x\sigma_z+\lambda y\sigma_x+\xi\sigma_y.
$$

At `xi=0`, the two surfaces form an exact conical intersection at the origin.  A
nonzero complex `sigma_y` SOC term opens a gap `2|xi|`.

## Complete Kramers-doublet model

With orbital Pauli matrices `tau` and spin Pauli matrices `sigma`,

$$
\mathbf V_D=
V_0\mathbf 1_4+\kappa x\,\tau_z\otimes\mathbf1
+\lambda y\,\tau_x\otimes\mathbf1
+\xi\,\tau_y\otimes\sigma_z.
$$

The antiunitary operator is

$$
\Theta=\mathbf1_{\mathrm{orbital}}\otimes(i\sigma_y)K,
\qquad \Theta^2=-1.
$$

The SOC term is time-reversal invariant and both eigenvalue pairs remain exactly
Kramers degenerate.  The validation residual is zero to displayed precision and the
maximum splitting is below `9e-19` hartree.

## Complete singlet/triplet model

The electronic basis contains two singlets followed by all three triplet
projections.  The singlet block carries the CI.  A fixed complex `2 x 3` SOC block
couples it to the triplet block, and Hermiticity fixes the reverse block.  Projector
ranks are exactly two and three and resolve the complete five-state identity.

This is an analytic topology and spin-completeness benchmark.  Its constant SOC
matrix is not claimed to be a quantitatively accurate molecular SOC surface.

## Numerical evidence

- Maximum norm drift on the primary trajectory: `0.0` at stored precision.
- Forward/backward phase-aligned error: approximately `4.50e-15`.
- Maximum edge-strip probability: approximately `8.23e-25`.
- Observed temporal orders: approximately `2.07` and `2.32`.
- Energy drift decreases by approximately four under each factor-two timestep
  refinement.
- Complete-doublet and complete-singlet/triplet smoke trajectories conserve norm to
  roughly `1e-15`.

## Boundary

No absorber or complex absorbing potential is included.  A calculation that allows
appreciable density to reach the periodic edge is not admitted as an exact reference.
