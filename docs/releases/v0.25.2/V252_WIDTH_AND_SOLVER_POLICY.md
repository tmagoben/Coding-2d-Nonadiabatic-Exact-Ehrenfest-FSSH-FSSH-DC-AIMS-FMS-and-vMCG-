# v0.25.2 width and solver policy

## Frozen defaults

| Quantity | Default |
|---|---:|
| Relative singular-value cutoff | `1e-10` |
| Absolute singular-value cutoff | `1e-12` |
| Maximum retained condition number | `1e10` |
| Compatible-null RHS tolerance | `2e-9` |
| Linear metric residual tolerance | `2e-9` |
| Nonlinear midpoint residual tolerance | `2e-10` |
| HYBR `xtol` | `1e-10` |
| Maximum function evaluations | `800` |
| Structural receipt tolerance | `2e-10` |
| Maximum per-step norm drift | `3e-8` |
| Minimum width | `1e-8` |
| Maximum width | `1e8` |
| Maximum absolute chirp | `1e8` |
| Maximum one-step `|Delta log(alpha)|` | `0.5` |

Log-width coordinates guarantee positivity mathematically. The explicit width and
chirp bounds prevent numerically meaningless yet positive states. The per-step log
gate limits a width ratio to `exp(0.5)` in either direction and detects a nonlinear
root that jumps to a remote branch.

## Nonlinear acceptance

The explicit tangent predictor initializes `scipy.optimize.root(method="hybr")`.
An endpoint is accepted only when HYBR reports success and the independently
recomputed full residual norm is below tolerance. The `xtol=1e-10` setting was chosen
after a constant-gauge solve reached machine residual at `1e-11` but HYBR reported
roundoff-scale stagnation; this is recorded as a solver-policy correction, not hidden.

## Metric acceptance

Full SVD exposes every retained and discarded direction. Admission requires a real
symmetric positive-semidefinite metric, at least one retained singular direction,
bounded retained condition number, compatible RHS projection onto discarded left
singular vectors, and a bounded reconstructed linear residual. No diagonal shift or
unreported regularizer is applied.

## Receipt trust chain

Validation binds model fingerprint, exact settings, signed time, packet count,
log-width midpoint, endpoint width/chirp domain, metric, RHS, velocity, SVD spectrum,
rank/nullity, nonlinear solver status and counts, residual vector/norm, norm, energy,
and step-chain continuity.

