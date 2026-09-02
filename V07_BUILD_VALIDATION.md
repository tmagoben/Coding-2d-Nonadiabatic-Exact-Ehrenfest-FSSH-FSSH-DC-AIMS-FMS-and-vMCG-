# v0.7 Build Validation Report

Validated on 2026-08-12.

## Compilation

All Python source files in the package, examples, and tests compiled successfully.

## Automated regression suite

```text
........................................................................ [ 91%]
.......                                                                  [100%]
79 passed in 3.07s
```

The cumulative suite includes all v0.1-v0.6 regression tests plus v0.7 tests for:

- flat and frustrated electronic gauge graphs;
- unitary polar links;
- gauge-invariant Wilson-loop spectra and traces;
- a discrete Berry phase of $-1$ around the v0.4 conical intersection;
- spanning-tree gauge transport;
- multi-start graph synchronization without removing physical holonomy;
- derivative-Hamiltonian reconstruction;
- arbitrary local $U(2)$ gauge covariance of electronic pair factors;
- gauge invariance of graph-Gaussian $S$ and $H$ matrices;
- generalized Cayley norm conservation;
- PySCF snapshot-graph construction using an injected many-electron overlap engine;
- deterministic TBF-center/pair-centroid graph topology.

## Conical-intersection graph regression

```text
CI lower-state gauge graph
--------------------------
number of nodes: 80
Wilson loop: (-1+0j)
Berry phase / pi: 1.0
tree-gauge objective: 3.9969161449628916
synchronized objective: 0.3477387490670639

The nontrivial W=-1 holonomy remains after gauge smoothing.
```

The synchronization objective decreases substantially while the Wilson loop remains
$-1$.  This is the intended behavior: finite-step/local gauge variation can be
redistributed, but the physical loop holonomy is not erased.

## Gauge-covariant graph-Gaussian regression

```text
Gauge-covariant graph-Gaussian matrix elements
-----------------------------------------------
max |S-S'|: 5.551115123125783e-16
max |H-H'|: 8.355534721610418e-17
Hermiticity residual H: 2.7755575615628914e-17
```

The matrices are invariant to machine precision under independently chosen random
$U(2)$ gauges at every electronic graph node.

## PySCF runtime status

PySCF is not installed in the build environment used for this release.

Therefore the real PySCF graph example

```bash
python examples/21_pyscf_tbf_centroid_graph.py
```

is included but was not numerically executed here.  Its graph edges call the v0.6
many-electron CASSCF overlap engine, whose algebra and API contract remain covered by
the inherited test suite.

## Scientific interpretation

Passing these tests establishes internal gauge covariance and numerical consistency of
the graph layer.  It does not imply that a real molecular graph is sufficiently dense,
that the selected electronic manifold is complete, or that nontrivial Wilson loops
should vanish.  Those are physical/convergence questions to be checked for each
application.
