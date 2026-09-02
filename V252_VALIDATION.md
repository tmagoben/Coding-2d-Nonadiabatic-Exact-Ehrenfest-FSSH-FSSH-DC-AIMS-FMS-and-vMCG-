# v0.25.2 validation

v0.25.2 adds 70 deterministic numerical gates and 25 adversarial/core gates. With
535 inherited v0.25.1 checks, the cumulative campaign contains 630 native Booleans.

## Canonical observations

| Evidence | Result |
|---|---:|
| Even maximum norm drift | `5.11e-15` |
| Odd maximum norm drift | `6.88e-15` |
| Even maximum energy drift (hartree) | `9.88e-14` |
| Odd maximum energy drift (hartree) | `4.54e-14` |
| Even maximum width change | `1.51e-4` |
| Odd maximum chirp change | `1.45e-3` |
| Largest odd reversal width/chirp error | `1.47e-13` |
| Permutation coefficient error | `1.39e-17` |
| Constant-gauge coefficient error | `2.00e-16` |
| Duplicate-packet rank/nullity | `6 / 6` of 12 |
| Duplicate-packet relative null forcing | `3.36e-16` |
| Largest canonical nonlinear residual | `<3.97e-16` |
| Harmonic endpoint width error at 20 au | `7.34e-9` |
| Harmonic endpoint chirp error at 20 au | `3.51e-8` |
| Refinement ratios | `0.25000008`, `0.24999887` |

The analytic overlap/Hamiltonian and complete adaptive metric/RHS are independently
reconstructed by dense-grid quadrature and numerical parameter tangents.

## Fail-closed controls

Controls reject disabled adaptive widths, spawning, pruning, coordinate-dependent
frames, multidimensional/full widths, real molecular-provider requests, non-SVD
solves, non-midpoint integration, wrong width coordinates, inverted width bounds,
indefinite metrics, incompatible null forcing, static providers, nonlinear failure,
width/chirp domain violations, excessive one-step log-width change, and tampered
endpoint, metric, RHS, velocity, singular spectrum, or nonlinear residual.

## Interpretation

This validates the mathematics and implementation for the exact released model
class. It is not a molecular accuracy benchmark or a basis-completeness argument for
general anharmonic dynamics.

