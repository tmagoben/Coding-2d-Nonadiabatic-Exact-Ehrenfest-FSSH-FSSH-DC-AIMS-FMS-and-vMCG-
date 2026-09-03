# v0.14 Release Notes

v0.14 is the **fully time-adaptive TDSE-defect control and complexity-audit release**.

The version intentionally skips a standalone v0.14 package: the planned v0.14
time-adaptive controller and the planned v0.14 complexity/budget work were integrated
into one release so the adaptation logic could be validated together with its
computational cost.

## New modules

```text
complexity_v14.py
fast_lvc_matrices_v14.py
residual_pruning_v14.py
defect_candidates_v14.py
adaptive_defect_dynamics_v14.py
v14_benchmark.py
```

## New numerical capabilities

```text
periodic TDSE-defect checks
hysteretic add/remove thresholds
adaptation cooldown
energy-conserving local dynamic candidate dictionaries
vectorized residual candidate ranking
zero-amplitude defect-triggered growth
exact leave-one-out nonorthogonal pruning loss
basis-budget replacement logic
emergency conditioning control
Hermitian half-build of exact S/H matrices
algorithmic timing/counter ledger
```

## Reference result

```text
initial basis: 10
final basis: 11
average basis: 10.925

projected-state dynamics error:
9.527804623132635e-05

target density error:
0.03330494031479218

population error:
0.028084897912098693

coherence phase error:
0.0028906431794148953

norm drift:
2.115581487549534e-06

maximum condition:
1470.7558920505405
```

## Adaptive event

```text
step: 10
time: 0.05
candidate: parent=9;target=0;pos=target_force:0.06;mom=nac;width=1.35

relative defect:
0.03238019095782191
    ->
0.029987957971150208

predicted capture fraction:
0.15067690241980206
```

## Complexity

The adaptive run performed:

```text
7931 Hermitian pair evaluations
versus
14531 ordered-pair equivalents
```

for a reduction of

```text
45.420 %
```

The observed runtime is dominated by exact Gaussian S/H and moving-basis T pair
algebra, not by the tiny 22-dimensional Cayley solve.

See `V14_ALGORITHM_COMPLEXITY.md`.

## Pruning stress test

```text
fractional loss: 0.0
condition before: 379319.12346481933
condition after: 67.30166373596352
improvement factor: 5636.103216600353
```

## Validation

All configured release acceptance checks pass.

The final cumulative automated suite reports:

```text
167 passed
```

Full details are recorded in `V14_BUILD_VALIDATION.md`.
