# v0.18 Convergence-Completeness Campaign

## Goal

The campaign is designed to distinguish four numerical questions that are often
collapsed into one:

```text
1. timestep discretization
2. Gaussian basis completeness
3. sparse graph truncation
4. adaptive growth-trigger sensitivity
```

Each coordinate is run in a **fresh Python process**.

This avoids process-history effects from previous FFT/sparse-solver work and ensures
that convergence points do not inherit mutable caches from one another.

## Exact references

Two exact-grid trajectories are propagated:

```text
original target packet
finite-Gaussian projected initial packet
```

Both use:

```text
exact grid dt = 0.0025
stored interval = 0.10
final time = 0.60
```

The conserved overlap between these two exact trajectories is used to separate initial
representation error from Gaussian propagation error.

## Canonical Gaussian coordinate

```text
dt = 0.005
max basis = 13
local omitted-score budget = 0.010
enrichment relative threshold = 0.015

defect control interval = 0.05 time units
minimum adaptation separation = 0.05
prune age = 0.10
cost horizon = 0.05
sampled audit interval = 0.10

candidate batch size = 16
```

All cadence quantities are converted to step counts only after `dt` is chosen.

## Basis axis

```text
Nmax = 10, 11, 12, 13
```

All other canonical controls are fixed.

## Timestep axis

```text
dt = 0.010, 0.005, 0.0025
```

The adaptive event times are held fixed in physical time.

Self-convergence is computed from the complex Gaussian wavefunctions, not from scalar
observables.

## Sparse-edge-budget axis

```text
B_local = 0.030, 0.010, 0.000
```

`B_local = 0` restores every locally scored omitted edge.

## Growth-trigger axis

```text
eta_enrich = 0.050, 0.035, 0.030, 0.025, 0.015
```

This deliberately spans:

```text
no enrichment
partial/delayed enrichment
full 13-Gaussian enrichment
plateau
```

## Audits

Normal trajectory checkpoints use sampled omitted-edge S/H/T scoring.

Full dense S/H/Snuc audits are retained only as initial and final sentinels.

## Machine-readable results

```text
results/v018_convergence_complete_campaign.json
```

The independent coordinate outputs are retained under:

```text
results/v018_partials/
```

This allows individual convergence points to be inspected without rerunning the full
campaign.
