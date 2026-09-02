# v0.17 Release Notes

v0.17 is the **sparse error-control release**.

## New modules

```text
edge_importance_v17.py
sparse_error_complexity_v17.py
error_controlled_sparse_dynamics_v17.py
sparse_error_budget_v17.py
v17_benchmark.py
```

## Main changes

```text
S/H/T edge importance replaces overlap-only final edge selection
score hysteresis
global local-score L2 budget
budget-driven edge promotion
periodic dense online S/H/Snuc audits
one-sided automatic threshold relaxation
edge-score convergence sweep
local-budget convergence sweep
actual graph-construction scaling benchmark
```

## Release result

```text
projected-state dynamics error:
0.00013361460054442858

target density error:
0.03333954068459557

population error:
0.02819941365898425

coherence phase error:
0.0029095064228609707

norm drift:
2.0053154399235495e-06
```

## Online controller

The deliberately aggressive initial graph fails its first matrix audit and
automatically relaxes once.

Final audit:

```text
S error:
0.005191742661052565

H error:
0.003962632349871911
```

Both satisfy the `0.006` release budget.

## Convergence

The final snapshot shows monotone error reduction under:

```text
decreasing edge-score threshold
decreasing local omitted-score budget
```

At the finest threshold:

```text
S error:
0.00026621461118036714

H error:
0.00015704048353805722
```

## Scaling

At `N=160` on the bounded-locality chain:

```text
pair-factorization reduction:
93.87 %

dense/sparse assembly speedup:
130.54 x
```

The exact S/H/T score-check exponent is approximately

```text
1.0557080719105292
```

while dense pair count is approximately

```text
1.97980963642741.
```

These are benchmark-specific scaling diagnostics.

## Scope

v0.17 remains an analytic-LVC sparse Gaussian reference implementation.

SOC is intentionally deferred until sparse/basis/dynamics convergence details are
better established.


## Automated regression

```text
200 passed
```

See `V17_BUILD_VALIDATION.md` for the complete validation record.
