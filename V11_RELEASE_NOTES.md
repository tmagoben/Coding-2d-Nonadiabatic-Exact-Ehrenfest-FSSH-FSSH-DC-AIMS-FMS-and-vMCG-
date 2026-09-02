# v0.11 Release Notes

Version 0.11 is the **basis-completeness and optimized-spawning-inspired release**.

## New modules

```text
gaussian_general.py
spa_matrix_elements_v11.py
moving_graph_gaussian_v11.py
optimized_spawning.py
managed_graph_aims_v11.py
basis_completeness.py
v11_benchmark.py
```

## New mathematical capabilities

- exact unequal-width Gaussian overlaps;
- exact unequal-width cross centroids and covariances;
- unequal-width kinetic matrix elements;
- unequal-width moving-basis matrix elements;
- unequal-width SPA0/SPA1 saddle points;
- energy-constrained child placement at shifted positions;
- local first-order coupling ranking;
- novelty-aware candidate selection;
- width-diverse children;
- multiple nonredundant children per event;
- descendant spawning with explicit lineage;
- canonical-basis participation diagnostics.

## Backward compatibility

The v0.10 managed propagator is retained unchanged as a public API.

v0.11 introduces

```python
run_basis_complete_graph_aims(...)
```

as a separate propagator.

This means all old benchmark results and unit tests remain reproducible rather than
silently changing under the same function name.

## Terminology

The release deliberately uses

> optimal-spawning-inspired local search

rather than

> implementation of optimal spawning.

The original Yang-Coe-Kaduk-Martinez method is a continuous nonlinear constrained
optimization designed to maximize parent-child Hamiltonian coupling at equal
classical energy.

v0.11 uses a finite local candidate search and a first-order coupling proxy.

## Main scientific question

Does a more diverse, better placed, multi-generation Gaussian basis repair the
basis-completeness failure diagnosed by v0.10?

The answer is recorded quantitatively in `V11_BENCHMARK_RESULTS.md`.


## Actual release benchmark

For the strong-CI benchmark:

```text
exact populations: [0.22600611046735578, 0.7739938895326441]
v0.11 populations: [0.23511178903234, 0.7648882109676599]

population error: 0.0128773741212
density error: 0.15991833275
exact purity: 0.676208196977
v0.11 purity: 0.643115591123
```

Relative to the v0.10 baseline:

```text
population-error improvement: 79.541 x
full-density-error improvement: 6.437 x
purity-error improvement: 9.025 x
```

The release passes the population and purity criteria but fails the stricter
full-density criterion. This points to electronic coherence/phase accuracy as the
next major target after the present basis-completeness improvement.


## Automated validation

```text
125 passed
```
