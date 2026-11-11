# v0.30.1 development validation

This report covers the recovered curved-frame modules and the new covariant
kinetic-plus-electronic action. It does not certify a new release.

## Current checks

- The saved v0.30.0 archive matches its recorded SHA-256 checksum.
- Recovered module suite before the new work: **36 passed**.
- Combined moving-frame, curved-basis, pre-ML, molecular-bridge, new Hamiltonian,
  and top-level v0.29.4/v0.29.5 integration suite: **131 passed** in 190.82 seconds.
- After tightening uniform-grid validation and adding two very-small-spacing
  regressions, the affected grid tests passed: **13 passed, 14 deselected**.
- The reproducible analytic-SOC example passed **8/8 numerical gates**.
- New Python modules compile, and `git diff --check` passes.

Runtime: Python 3.12.14, NumPy 2.3.5, SciPy 1.17.0, with BLAS/OpenMP threads set
to one. The pytest run emits a deprecation warning from the installed AnyIO pytest
plugin; no project test is skipped by these focused commands.

## Measured evidence

| Quantity | Result |
|---|---:|
| Full-Hamiltonian moving-gauge covariance error | 1.1075e-16 |
| Hermiticity residual | 5.9797e-18 |
| Maximum NAC/overlap link disagreement | 6.1022e-16 |
| Finite-grid propagation gauge error | 2.1048e-16 |
| Finite-grid norm drift | 4.4409e-16 |
| Finite-grid energy drift, hartree | 2.1684e-19 |
| Observed kinetic convergence orders | 1.9710, 1.9927 |
| TDVP metric rank | 12 |
| Retained metric condition number | 96.7636 |
| Relative SVD equation residual | 1.1976e-15 |

The continuum kinetic comparison uses grid spacings 0.4, 0.2, and 0.1 bohr.
Relative errors are 0.0438311, 0.0111802, and 0.00280923. The finite-grid exponential
is an independent reference calculation; it is not a Gaussian time trajectory.

Machine-readable results are in `results/v0301_molecular_hamiltonian_evidence.json`.
Source recovery hashes are in `results/v0301_recovery_manifest.json`.

## Cumulative baseline

The untouched GitHub baseline is commit `7e7f53f0ec63d496e1e7ec6bf605e557b22b1997`.
Its existing [October 10 cumulative CI run](https://github.com/tmagoben/Coding-2d-Nonadiabatic-Exact-Ehrenfest-FSSH-FSSH-DC-AIMS-FMS-and-vMCG-/actions/runs/38013590799)
finished with **22 failed, 617 passed, 16 skipped**. Its failures include missing
PySCF 2.13.1, absent historical result JSON files, and two machine-precision
assertions. The development-integration job on that same commit passed.

A local cumulative run was started before recovery or implementation changes.
Its result will be recorded separately when complete. The focused passing results
above must not be presented as a green cumulative or release campaign.

## Acceptance boundary

Validated here: the declared constant-diagonal-mass projected operator, discrete
Hermiticity/covariance, analytic kinetic refinement, and its connection to the
full-parameter TDVP snapshot. Outstanding: implicit-midpoint Gaussian trajectories,
adaptive-event integration, dataset schema reconciliation and campaign certification,
and live molecular-SOC validation.
