# v0.2 Validation Contract

Version 0.2 adds genuinely new nonadiabatic structure, so its acceptance criteria are
stricter than "the spawning example transfers population."

- The exact two-state split operator must conserve total norm.
- Adiabatic eigenvectors must be sign aligned on the 1D path.
- $d_{01}=-d_{10}$.
- The covariant adiabatic Hamiltonian
  $-\frac{1}{2M}(\partial+d)^2+E$ must reproduce the diabatic Hamiltonian after basis
  transformation within finite-difference error.
- Adiabatic Gaussian $S$ and $H$ matrices must be Hermitian to stated quadrature/
  derivative tolerance.
- An energy-conserving child must satisfy the local classical energy equation.
- Adding a zero-amplitude child must leave the represented wavefunction unchanged.
- Spawning is called a foundation/prototype unless full FMS/AIMS spawning optimization
  and convergence machinery is implemented.
