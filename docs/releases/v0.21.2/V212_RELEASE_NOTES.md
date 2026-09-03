# v0.21.2 Release Notes

v0.21.2 is the **pre-SOC integration-hardening release**.

It deliberately contains **no physical spin-orbit Hamiltonian**.

## Added

- unequal-width complex block Gaussian $S/H/T$ using the exact general Gaussian algebra;
- representation-neutral local block mean-field nuclear guidance;
- self-consistent velocity-Verlet + coefficient predictor/corrector propagation;
- optional adaptive step-boundary block-basis actions;
- exact zero-electronic-block Gaussian insertion;
- metric-projected block pruning using the Schur complement;
- generic Hermitian electronic observable matrices and expectations;
- full-subspace continuity diagnostics integrated with an indexed provider wrapper;
- pre-SOC complex-dtype static/runtime audit;
- explicit package discovery and a clean editable-install regression guard;
- v0.21.2 acceptance campaign and regression tests.

## Canonical hardening results

```text
unequal-width S covariance error: 2.311047820503235e-16
unequal-width H covariance error: 3.1082368759020804e-16
unequal-width T covariance error: 1.6051159987815854e-16
```

Self-consistent gauge-equivalence convergence:

```text
orders: [2.002130949292044, 1.995625423493078]
finest coefficient error: 1.2184279933242558e-13
finest maximum norm drift: 4.085620730620576e-14
```

Generic observable gauge error:

```text
1.734723475976807e-18
```

Complex dtype audit:

```text
passed = True
suspicious casts = []
```

## Scientific scope

Use:

> pre-SOC-hardened, complex representation-neutral block-sparse Gaussian nonadiabatic
dynamics framework with unequal-width pair algebra, coefficient-coupled nuclear
guidance, block adaptive lifecycle, subspace continuity diagnostics, and generic
electronic observables.

Do not describe v0.21.2 as production AIMS or as SOC dynamics.
