# v0.4 Validation Contract

Version 0.4 is accepted only if its topology, Gaussian algebra, exact reference, and
coupled spawned-basis propagation all pass independent checks.

## Conical-intersection model

- $V_d(x,y)$ is Hermitian.
- $E_+(0,0)=E_-(0,0)$.
- Analytic adiabatic energies match `numpy.linalg.eigvalsh`.
- Analytic gradients match centered finite differences away from the CI.
- Vector NACs satisfy $\mathbf d_{+-}=-\mathbf d_{-+}$.

## Geometric phase / gauge

- Numerical line integration of $\mathbf d_{-+}$ around a closed loop gives $\pi$
  for one winding.
- Continuous sign transport of a real adiabatic eigenvector around the CI returns a
  final overlap near $-1$ with the initial vector.
- Unitary Procrustes transport correctly removes an arbitrary rotation inside a
  selected subspace.

## Multidimensional Gaussians

- Grid quadrature gives unit norm.
- Coordinate covariance agrees with $\frac12A^{-1}$.
- Analytic equal-width overlap agrees with 2D quadrature.
- Analytic kinetic action agrees with spectral differentiation on a localized packet.
- Heller width matrices remain symmetric.
- For a 2D harmonic oscillator, the coherent-state width remains constant and the
  center follows the analytic classical solution.

## Exact 2D reference

- Two-state FFT propagation conserves total norm.
- Electronic populations sum to one.
- Reducing $\Delta t$ produces convergent results.

## Coupled spawned Gaussian basis

- $S$ and $H$ are Hermitian within quadrature/finite-difference tolerance.
- $\dot S=T+T^\dagger$ for the moving basis.
- A zero-amplitude spawned child leaves the wavefunction unchanged at insertion.
- A real NAC-rescaled child conserves local trajectory energy.
- Coefficient propagation conserves $C^\dagger S C$ to the expected integration
  tolerance before/after basis growth.
- A fixed set of inputs gives deterministic spawn decisions.
- Overlap blocking prevents repeated insertion of a nearly redundant target TBF.

Passing these tests establishes internal consistency. It does not by itself establish
chemical convergence for a real molecule.
