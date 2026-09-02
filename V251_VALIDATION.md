# v0.25.1 validation

The release defines 55 deterministic numerical gates and 20 core/adversarial gates.
Together with all 460 v0.25.0 gates, the cumulative contract contains 535 native
Boolean checks.

## Positive evidence

| Evidence | Canonical result |
|---|---:|
| Even-state maximum norm drift | `7.24e-14` |
| Odd-state maximum norm drift | `3.06e-14` |
| Even maximum energy drift (hartree) | `6.97e-13` |
| Odd maximum energy drift (hartree) | `2.72e-13` |
| Largest even/odd reversal coefficient error | `6.94e-17` |
| Gaussian-permutation coefficient error | `5.72e-17` |
| Constant-gauge coefficient error | `6.21e-17` |
| Duplicate-packet metric rank/nullity | `4 / 4` of 8 |
| Duplicate-packet relative null forcing | `5.74e-17` |
| Duplicate-packet relative solve residual | `1.17e-15` |
| Largest canonical nonlinear residual | `<4.91e-16` |
| Refinement ratios | `0.24999991`, `0.25000056` |

The production analytic overlap/Hamiltonian and full tangent metric/RHS are also
compared against independent dense-grid quadrature in the core tests.

## Negative controls

The cumulative release fails closed on adaptive widths, spawning, pruning,
coordinate-dependent frames, multidimensional requests, real molecular-provider
admission, non-SVD metric solves, non-midpoint integration, incompatible null-space
forcing, indefinite metrics, static-provider intake, nonlinear nonconvergence, and
tampered endpoints, metric, RHS, velocity, singular spectrum, or residual.

## Interpretation

These results establish the algebra, solver, invariances, and receipt integrity of
the released model class. They are not evidence for accuracy on a general molecule,
for adaptive-basis completeness, or for long-time production performance.

