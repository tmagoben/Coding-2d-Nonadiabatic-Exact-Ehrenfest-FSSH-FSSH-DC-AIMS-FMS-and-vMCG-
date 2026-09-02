# v0.9 Release Notes

## Theme

**Convergence before complexity.**

v0.9 does not add another representation of nonadiabatic dynamics.  It adds the tools
needed to decide whether the existing graph-Gaussian approximation is numerically
under control.

## New modules

```text
gaussian_dynamics/spa_matrix_elements.py
gaussian_dynamics/basis_management.py
gaussian_dynamics/adaptive_spawning.py
gaussian_dynamics/convergence.py
gaussian_dynamics/managed_graph_aims.py
gaussian_dynamics/exact_benchmark.py
gaussian_dynamics/benchmark_suite.py
```

## New examples

```text
26_spa0_vs_spa1.py
27_basis_pruning.py
28_adaptive_spawning.py
29_managed_vs_exact_ci.py
30_timestep_refinement.py
```

## Main scientific additions

1. explicit SPA0 and first-order electronic Taylor matrix elements;
2. exact first-moment Gaussian correction through the complex pair centroid;
3. overlap-eigenvalue conditioning diagnostics;
4. Hilbert-space projection when pruning redundant TBFs;
5. cumulative pruning-loss budget;
6. timestep-aware integrated nonadiabatic coupling exposure for spawning;
7. managed graph-AIMS propagation with SPA order and basis conditioning controls;
8. exact 2D adiabatic-population benchmark;
9. reusable convergence/refinement utilities.

## Terminology

The `order=1` matrix-element implementation is deliberately described as an
**SPA1 electronic Taylor layer** and not as a complete production AIMS-SPA1
Hamiltonian.
