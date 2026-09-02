# Reference-First Multidimensional CI+SOC Dynamics

v0.26.0 follows a strict validation order:

1. define a Hermitian fixed-frame quadratic CI+SOC Hamiltonian;
2. converge an independent two-dimensional FFT-grid trajectory;
3. evaluate the same initial spinor with the analytic Gaussian ansatz;
4. propagate fixed-basis and controlled adaptive TDVP trajectories;
5. compare phase-aligned wavefunctions and reduced electronic densities;
6. admit a release claim only after reductions, symmetries, conditioning, and
   adversarial controls pass.

The exact and Gaussian solvers intentionally share only the Hamiltonian model and
initial physical state.  They do not share propagation algebra.

The release provides complete Kramers-doublet and singlet/triplet model spaces, but
these are analytic validation systems.  PySCF remains an electronic-structure source
elsewhere in the repository; v0.26.0 does not call PySCF inside the multidimensional
trajectory loop.

For derivations and numerical evidence, see the root-level v0.26.0 documents.
