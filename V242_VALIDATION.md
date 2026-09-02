# v0.24.2 validation

## Canonical real calculation

- PySCF 2.13.1, CPython 3.12, Linux x86-64.
- OH radical, neutral charge, doublet sector, STO-3G.
- ROHF followed by equal-weight three-root SA-CASSCF(5e,4o).
- Center bond length `1.83256418024373` bohr.
- Centered displacements `0.08`, `0.04`, and `0.02` bohr.
- Seven independently converged geometry snapshots.
- Three spin-free roots and six complete doublet microstates per snapshot.

## Numerical observations

The direct-JK and explicit center SOMF paths agree with maximum absolute error
`1.7763568394002505e-15`. Across all six center-to-endpoint overlaps, the minimum
retained singular value is approximately `0.9968178664` and the maximum is
approximately `0.9999117824`. Maximum polar-transport nonunitarity is below
`1.7e-15`; the audited matrix/symmetry residual is below `6.0e-12`.

The observed fine-to-coarse change ratios are approximately

| Component | Ratio | Ideal centered-difference ratio |
|---|---:|---:|
| `K_spin_free` | 0.252553 | 0.25 |
| `K_soc` | 0.257496 | 0.25 |
| `K_total` | 0.252553 | 0.25 |

The relative Richardson estimates are about `1.88e-3` for the spin-free/total
component and `3.33e-4` for SOC. The finest SOC derivative Frobenius norm is about
`1.378e-4` hartree/bohr, so the convergence signal is nonzero.

## Acceptance structure

- 315 gates inherited unchanged from v0.24.1.
- 60 real PySCF connected-geometry runtime gates.
- 25 core/adversarial gates.
- 85 new gates and 400 cumulative gates.

Negative controls reject spectral expansion, rank loss, ambiguous degenerate
root-phase assignment, non-Hermitian derivatives, direct-JK/oracle disagreement,
inflated capability claims, nondecreasing step ladders, endpoint-fingerprint
substitution, and equal-and-opposite component tampering.

## Interpretation

This validates implementation consistency and one-coordinate numerical convergence.
It does not validate a general molecular SOC model, basis/method convergence,
spectroscopic accuracy, full Cartesian derivatives, analytic derivatives, continuous
physical SOC-spinor derivative connections, real mixed multiplicity, or trajectory
admission.
