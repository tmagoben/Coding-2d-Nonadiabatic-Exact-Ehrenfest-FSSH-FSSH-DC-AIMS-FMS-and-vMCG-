# v0.7 Validation Contract

## Graph connection

- polar links are unitary;
- reverse links equal the Hermitian conjugate;
- coefficient transport composes correctly along paths;
- a flat synthetic $U(2)$ graph can be globally trivialized.

## Wilson loops

- Wilson-loop eigenvalues and trace are invariant under arbitrary local $U(2)$ gauge
  transformations;
- a synthetic $U(1)$ triangle with phase $\pi$ retains $W=-1$ under synchronization;
- the lower-state ring around the analytic v0.4 conical intersection gives $W=-1$.

## Gauge synchronization

- spanning-tree transport makes all tree links equal to identity;
- the synchronization objective does not increase relative to the tree gauge in the
  frustrated $U(1)$ regression;
- nontrivial holonomy remains nontrivial after synchronization.

## Electronic operators

- derivative-Hamiltonian matrices reconstructed from energies, gradients, and NACs are
  Hermitian for real adiabatic input;
- pair electronic overlap, potential, and derivative-Hamiltonian factors are invariant
  under independent random local $U(2)$ gauges.

## Graph-Gaussian layer

- static graph-Gaussian overlap and Hamiltonian matrices are Hermitian for a symmetric
  pair-reference choice;
- the full matrices are unchanged by arbitrary independent node gauges;
- generalized Cayley propagation preserves $C^\dagger SC$.

## PySCF graph helper

- snapshot graph construction accepts an injected many-electron overlap engine;
- edge singular-value/unitarity diagnostics are produced;
- TBF-center/pair-centroid connectivity is generated deterministically.

## Full regression

All v0.1-v0.6 tests must remain passing.

## Research-use requirement

For a real PySCF graph, additionally check convergence with respect to

1. graph-node spacing;
2. number of electronic states in the overlap manifold;
3. active space;
4. edge selection;
5. pair-centroid placement;
6. synchronization versus unsynchronized observables;
7. loop holonomy under graph refinement.
