# Full Correlated Width Matrices

v0.27.0 closes the main geometric limitation of v0.26.0. A coordinate-diagonal
width cannot remain diagonal under a general rotation and cannot describe a tilted
or correlated Gaussian ellipse. The new packet uses full symmetric matrices
$\Gamma$ and $B$, with $\Gamma=\exp(E)$.

The implementation uses an orthonormal symmetric packing (`svec`) so coordinate
rotations remain orthogonal maps in the nonlinear parameter vector. Analytic
multivariate moments through degree four then support the complete McLachlan metric
without numerical differentiation or grid quadrature.

The release is checked against two independent references:

1. direct dense-grid/FFT integration of correlated overlap and Hamiltonian elements;
2. the exact matrix Riccati equations for a Gaussian in a rotated quadratic
   Hamiltonian.

The same representation is used by controlled spawning, projection, pruning, merge,
and dormant activation. Nondegenerate width eigenvectors provide intrinsic candidate
directions. Degenerate directions fail closed because no unique physical principal
axis exists.

See the root-level v0.27.0 derivation, covariance, lifecycle, validation, complexity,
and architecture documents for the frozen equations and gates.
