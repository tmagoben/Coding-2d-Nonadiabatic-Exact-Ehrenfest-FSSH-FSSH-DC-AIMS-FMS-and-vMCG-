# v0.20 Release Notes

v0.20 is the **end-to-end sparse molecular machinery release**.

## New modules

```text
indexed_molecular_provider_v20.py
sparse_molecular_matrices_v20.py
sampled_molecular_audit_v20.py
sparse_molecular_dynamics_v20.py
v20_benchmark.py
```

## Main changes

```text
buffered KD-tree electronic tracking cache
uid-based persistent molecular locality graph
nuclear-overlap geometric pre-screen
molecular S/H/T edge scoring
global local omitted-score budget
active-edge sparse molecular S/H/T assembly
sparse metric-compatible moving-basis connection
sparse midpoint/Cayley propagation
sampled molecular omitted-edge audits
automatic search-floor relaxation
independent dense-sentinel electronic caches
molecular threshold and budget convergence campaigns
bounded-locality molecular scaling campaign
```

## Canonical result

```text
20 TBFs
190 possible pairs
36 active edges
average sparsity = 81.05 %

dense-metric coefficient error = 0.000670432070175427
norm drift = 8.881784197001252e-16
new electronic-point reduction = 49.30 %
```

## Final dense sentinel

```text
S error = 0.002954918852389199
H error = 0.0016775567668417452
T error = 0.01718013850955588
```

## Scaling

At `N=160`:

```text
active edges = 152
exact molecular pair checks = 317
dense pairs = 12720
pair-check reduction = 97.51 %
matrix sparsity = 98.81 %
```

Fitted molecular pair-check exponent:

```text
1.0324491201728092
```

Formal dense pair exponent:

```text
2.0213244854837162
```

## Online controller

A deliberately over-aggressive search floor misses an important molecular edge.

The sampled audit detects it and changes

```text
0.90 -> 0.45
```

before the immediate re-audit passes.

## PySCF status

PySCF remains unavailable in the build environment.

The sparse molecular machinery is connected to the existing raw SA-CASSCF snapshot and
many-electron overlap bridge, but runtime validation here uses the deterministic
Cartesian LVC molecular backend only.

## Release status

```text
passed = True
```

The cumulative automated regression suite reports:

```text
231 passed
```

Full details are recorded in `V20_BUILD_VALIDATION.md`.

## Scope

v0.20 finishes the sparse molecular software architecture.

It does not claim production AIMS matrix elements, a real PySCF runtime trajectory, or
spin-orbit-coupled dynamics.
