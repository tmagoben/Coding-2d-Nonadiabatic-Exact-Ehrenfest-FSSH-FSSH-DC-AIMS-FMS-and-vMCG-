# v0.14 Time-Adaptive Defect Campaign

## 1. Purpose

This campaign tests whether the instantaneous v0.13 TDSE-defect primitive can be
embedded in a stable time-dependent basis controller.

It is not tuned to minimize one final observable at any cost.

The controller must obey explicit:

```text
defect thresholds
basis-size budget
conditioning budget
candidate-capture threshold
adaptation cooldown
pruning-loss budget
```

## 2. Initial state

The intended exact initial state is projected into a 10-Gaussian spinor-complete bank
using the v0.13 residual/density-screened initial-basis algorithm.

The adaptive run and the exact projected-state reference begin from the same projected
wavefunction.

## 3. Propagation settings

```text
dt = 0.005
steps = 120
final time = 0.6
defect diagnostic grid = 40 x 40
defect interval = 10 steps
```

## 4. Controller settings

```text
add threshold = 0.020
remove threshold = 0.006
minimum candidate capture fraction = 0.003

minimum basis = 8
maximum basis = 11

adaptation separation = 10 steps
minimum prune age = 20 steps
pruning patience = 2 checks

ordinary prune loss budget = 5e-7
replacement prune loss budget = 5e-7

candidate position shifts = 0, +0.06, -0.06
candidate width scales = 0.75, 1.0, 1.35
candidate momentum directions = NAC, momentum
```

## 5. Growth

At each high-defect checkpoint:

1. generate local energy-conserving candidates around all current Gaussian centers;
2. include same- and other-guidance-surface candidates;
3. evaluate the candidate Gaussian grid values;
4. project candidates out of the current nuclear Gaussian span;
5. contract every candidate with the TDSE defect;
6. reject near-dependent/ill-conditioned candidates;
7. accept the largest residual-capture fraction above threshold;
8. insert with zero two-state coefficient;
9. rebuild the matrices;
10. recompute the defect to verify improvement.

## 6. Pruning

When the defect is persistently low, or when replacement/conditioning control requires
it:

1. compute one inverse/solve of the nuclear overlap matrix;
2. obtain every exact leave-one-out loss;
3. protect recently added basis functions;
4. remove only a Gaussian below the configured loss budget;
5. project the old state into the retained basis;
6. rebuild and renormalize;
7. record loss and conditioning.

The primary release trajectory does not need a real pruning event.

A separate deterministic pruning stress test validates the operation.

## 7. Complexity accounting

Every propagation run stores a `complexity` dictionary.

The key distinction is between:

```text
asymptotic complexity
```

and

```text
observed small-benchmark wall time.
```

At the current basis size, pairwise Gaussian algebra dominates runtime even though the
dense coefficient solve is asymptotically cubic.

See `V14_ALGORITHM_COMPLEXITY.md`.

## 8. Reproduction

Run:

```bash
python examples/58_recompute_v014_campaign.py
```

The canonical output is:

```text
results/v014_time_adaptive_defect_campaign.json
```
