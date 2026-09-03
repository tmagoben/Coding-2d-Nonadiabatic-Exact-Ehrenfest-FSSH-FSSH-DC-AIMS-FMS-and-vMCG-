# v0.16 Release Notes

v0.16 is the **persistent locality graph and sparse propagation release**.

## New modules

```text
locality_graph_v16.py
sparse_pair_matrices_v16.py
electronic_cost_v16.py
local_cost_aware_v16.py
sparse_complexity_v16.py
sparse_adaptive_dynamics_v16.py
v16_benchmark.py
```

## Main additions

```text
rigorous conservative overlap upper bound
safe KD-tree global-radius screening
persistent uid-keyed graph topology
enter/exit edge hysteresis
sparse CSR S/H/T matrices
sparse metric-compatible moving-basis connection
sparse Cayley solve
local-degree adaptation cost
electronic geometry-cache cost term
dense endpoint sparse-matrix audit
bounded-locality scaling benchmark
```

## Primary physical result

```text
projected-state dynamics error:
0.00013361460054812487

target density error:
0.03333954068459046

population error:
0.028199413658981914

coherence phase error:
0.0029095064228781115

norm drift:
2.0053154308197207e-06
```

## Sparse approximation

```text
average graph sparsity:
5.38 %

final rho difference vs v0.15:
0.0001631265880682423

final relative S error:
0.005191742661052565

final relative H error:
0.003962632349871911
```

## Scaling result

At $N=80$ on the bounded-locality Gaussian chain:

```text
edge fraction:
0.04968354430379747

pair-factorization reduction:
92.69 %

dense/sparse assembly speedup:
32.56 x
```

The fitted active-edge exponent is approximately

```text
1.0425836916313385
```

while dense canonical pair count has exponent approximately

```text
1.9737662900529327.
```

These exponents are benchmark-specific.

## Validation

All configured v0.16 release acceptance checks pass.

The final cumulative automated suite reports:

```text
189 passed
```

Full details are recorded in `V16_BUILD_VALIDATION.md`.
