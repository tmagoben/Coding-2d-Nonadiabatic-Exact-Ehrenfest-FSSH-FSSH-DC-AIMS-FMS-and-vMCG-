# v0.15 Release Notes

v0.15 is the **shared Gaussian-pair cache, incremental matrix, and cost-aware
adaptation release**.

## New modules

```text
pair_cache_v15.py
complexity_v15.py
cost_aware_adaptation_v15.py
defect_candidates_v15.py
adaptive_defect_dynamics_v15.py
v15_benchmark.py
```

## Main numerical changes

```text
one multi-RHS pair solve for overlap/centroid/covariance
canonical-pair conjugate reversal
cached S/H assembly
cached ordered moving-basis T assembly
endpoint cache reuse during TDSE-defect checks
candidate-cache reuse after acceptance
incremental S/H append
zero-integral S/H/cache pruning
cost-aware residual candidate utility
expanded cache/factorization ledger
```

## Reference physics

```text
projected-state dynamics error:
9.527804623556872e-05

target density error:
0.03330494031478426

population error:
0.028084897912094255

coherence phase error:
0.0028906431794135244

norm drift:
2.115581485107043e-06
```

Maximum v0.14 reference-metric difference:

```text
8.412825991399586e-12
```

## Performance

```text
v0.14 factorization-equivalent baseline:
103103

v0.15 propagation pair factorizations:
15675

factorization reduction:
84.797 %

saved timing speedup:
2.683 x
```

The wall-time number is diagnostic. The pair-factorization removal is the portable
algorithmic statement.

## Cost-aware event

```text
capture fraction:
0.1506769024198026

utility:
0.24210065440802597

normalized incremental cost:
0.6223729662698816

estimated horizon seconds:
0.03330978168687729
```

For this benchmark the cost-aware and residual-only selectors choose the same candidate.

## Validation

All configured v0.15 release acceptance criteria pass.

The final cumulative automated suite reports:

```text
178 passed
```

Full details are recorded in `V15_BUILD_VALIDATION.md`.
