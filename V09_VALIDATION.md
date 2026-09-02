# v0.9 Validation Contract

v0.9 adds convergence management, so its tests are organized by approximation layer.

## SPA layer

- SPA0 returns \(F(q_c)S_{ij}\).
- SPA1 exactly integrates a linear scalar field between equal-width Gaussians.
- Graph SPA1 matrices remain Hermitian with symmetric pair references.
- SPA0/SPA1 difference is reported rather than hidden.

## Basis management

- exact duplicate removal has zero projection loss;
- near-duplicate pruning decreases overlap condition number;
- pruning is rejected when the configured projection-loss budget would be exceeded;
- canonical orthogonalization satisfies \(X^\dagger SX=I\).

## Adaptive spawning

- integrated coupling exposure is dimensionless;
- constant-coupling trigger time is stable under timestep refinement;
- child insertion still starts with zero coefficient;
- child momentum remains energy conserving through the inherited v0.8 rule.

## Time propagation

- generalized norm remains stable in the managed graph-AIMS benchmark;
- the metric-compatible basis connection remains inherited from v0.8;
- spawn and prune operations are explicitly logged.

## Exact reference

- exact 2D split-operator norm remains one;
- final adiabatic populations sum to the exact norm;
- grid/timestep refinement remains available independently of the Gaussian calculation.

## Convergence harness

- observed-order utility reproduces a synthetic second-order sequence;
- managed timestep studies return final population vectors at every refinement level;
- exact-vs-managed comparison returns a population error, not merely two plots.

## Full regression

All v0.1-v0.8 tests must continue to pass.

## What a research benchmark should report

Before claiming convergence, record at minimum:

```text
exact grid shape / box / dt
Gaussian dt
SPA order
spawn action threshold
maximum basis size
number and times of spawns
number and loss of pruning events
max cond(S)
max norm drift
final population vector
population error vs exact reference
```

A single passing unit test is not a convergence study.
