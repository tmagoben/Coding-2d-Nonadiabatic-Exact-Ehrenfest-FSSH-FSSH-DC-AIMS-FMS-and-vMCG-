# v0.25.1 metric and nonlinear-solver policy

## Frozen defaults

| Quantity | Default |
|---|---:|
| Relative singular-value cutoff | `1e-10` |
| Absolute singular-value cutoff | `1e-12` |
| Maximum retained condition number | `1e10` |
| Compatible-null RHS tolerance | `2e-9` |
| Linear metric residual tolerance | `2e-9` |
| Nonlinear midpoint residual tolerance | `2e-10` |
| Nonlinear `xtol` | `1e-11` |
| Maximum function evaluations | `600` |
| Structural receipt tolerance | `2e-10` |
| Maximum per-step norm drift | `2e-8` |

The broad production thresholds are safety gates, not the achieved validation
accuracy. Canonical evidence is much tighter: full-rank condition number is about
`5.8e2`, duplicate-packet null forcing is about `5.7e-17`, accepted nonlinear
residuals are below `5e-16`, and short-trajectory norm drift is below `8e-14`.

## Why SVD rather than normal equations

The physical tangent metric is positive semidefinite, not guaranteed positive
definite. Forming normal equations or blindly regularizing would square conditioning
or silently invent dynamics in redundant directions. Full SVD exposes the retained
and discarded spaces directly, allows an explicit compatibility test, and returns
the minimum-norm velocity when the null forcing is physical.

Rank deficiency is not automatically an error. It is admitted only when the RHS has
negligible projection into the discarded left-singular space. An incompatible RHS,
an indefinite metric, no retained direction, excessive retained conditioning, or an
excessive linear residual all fail closed.

## Why implicit midpoint rather than Verlet

The unknowns are coupled noncanonical variational parameters, and their metric
depends on the current Gaussian state. Independent kick/drift updates do not solve
this system. Implicit midpoint advances the full TDVP vector field together, is
second order and self-adjoint, and admits a direct residual receipt. Velocity Verlet
remains valid for the narrower constant-mass canonical system in v0.25.0; it is not
used by the v0.25.1 multi-Gaussian solver.

## Receipt binding

Validation reconstructs the midpoint from the stored endpoints, rebuilds the exact
metric and RHS, recomputes the SVD solution and spectrum, reconstructs the nonlinear
residual, and recomputes endpoint norm and energy. It also binds packet count,
unchanged widths, model fingerprint, settings, signed time, and trajectory chain
continuity. Tampering with any of these fields is rejected.

