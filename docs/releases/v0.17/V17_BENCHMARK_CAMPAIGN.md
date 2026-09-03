# v0.17 Sparse Error-Control Benchmark Campaign

## Purpose

The v0.17 campaign asks:

> Can sparse Gaussian graph decisions be controlled by actual S/H/T importance and an
> explicit matrix-error budget rather than by overlap locality alone?

The campaign deliberately starts from a graph that is **too aggressive** and requires
the online controller to repair it.

## Primary benchmark

The physical problem is unchanged from v0.13-v0.16:

```text
2D analytic LVC conical-intersection model
initial residual-selected Gaussian basis = 10
mass = 5 a.u.
q0 = (-0.60, 0.25)
p0 = (10.0, 0.0)
final time = 0.60
Gaussian propagation dt = 0.005
exact-grid reference dt = 0.0025
```

The independent TDSE-defect grid remains `40 x 40`.

## Initial graph policy

```text
enter score = 0.060
exit score = 0.030
geometric search overlap floor = 1e-5

overlap weight = 1.0
Hamiltonian weight = 0.20
time-connection weight = 1.0

local omitted-score L2 budget = 0.08
```

These score thresholds are intentionally aggressive.

## Dense online audit

Every 20 steps, plus the initial checkpoint, v0.17 performs a complete dense pair
rebuild and compares:

```text
S
H
Snuc
```

against the current sparse matrices.

The accepted budget is

```text
relative S error <= 0.006
relative H error <= 0.006
relative Snuc error <= 0.006
```

On failure, all graph thresholds are relaxed by a factor `0.5`.

Online tightening is not permitted.

## Adaptive Gaussian basis

The inherited TDSE-defect basis-growth logic remains active.

The basis can grow from 10 to at most 11 TBFs.

The edge-error controller and TBF residual controller are separate:

```text
edge controller:
    decides which pair blocks exist

basis controller:
    decides which Gaussian functions exist
```

## Final snapshot convergence sweeps

Two separate dense-audited sweeps are performed.

### Edge-score sweep

```text
0.12
0.08
0.06
0.04
0.03
0.02
0.01
```

The global local-score budget is disabled so this sweep isolates the per-edge score.

### Local-score-budget sweep

At fixed nominal enter score `0.06`:

```text
unbounded
0.10
0.08
0.05
0.03
0.01
0.00
```

This isolates the effect of the global omitted-score proxy.

## Construction-scaling benchmark

The actual v0.17 S/H/T scoring graph is tested on bounded-locality Gaussian chains at

```text
N = 20, 40, 80, 160
```

The campaign records:

```text
KD-tree spatial pairs
pair-specific bound screens
exact S/H/T score checks
active graph edges
pair factorizations
dense canonical pair count
sparse H density
dense and sparse assembly time
```

Log-log exponents are descriptive diagnostics for this benchmark only.

## Reproduction

Run:

```bash
python examples/77_recompute_v017_campaign.py
```

Canonical output:

```text
results/v017_sparse_error_control_campaign.json
```
