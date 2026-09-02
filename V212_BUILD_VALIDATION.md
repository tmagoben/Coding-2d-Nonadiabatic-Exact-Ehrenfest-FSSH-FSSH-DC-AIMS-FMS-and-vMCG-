# v0.21.2 Build Validation

Validated on 2026-08-13.

## Source validation

```text
371 Python files parsed successfully with Python AST.
```

No unresolved release-template placeholders remain.

## Clean-install validation

The documented isolated editable-install path was rerun from the packaged source:

```text
python -m venv .venv
.venv/bin/python -m pip install -e ".[dev]"
```

The install succeeds. Package discovery is explicit in `pyproject.toml`, so the
top-level numerical `results/` directory is retained as release data without being
misclassified as a second Python package. A dedicated packaging regression test guards
this configuration.

## Automated regression suite

Pytest collected:

```text
246 tests
```

The final suite was verified in two deterministic groups to avoid repeatedly executing
the expensive historical release campaigns inside one monolithic run:

```text
242 non-release regression tests passed
4 versioned release-acceptance tests passed
-------------------------------------------
246 / 246 tests passed
```

The four release tests are:

```text
test_v019_release.py
test_v020_release.py
test_v021_release.py
test_v0212_release.py
```

The v0.21.2-specific focused suite also reports:

```text
9 passed
```

## Canonical campaign

```text
results/v0212_pre_soc_hardening_campaign.json
```

Release acceptance:

```text
16 / 16 configured checks passed
passed = True
```

## Unequal-width complex block algebra

```text
S covariance error = 2.311047820503235e-16
H covariance error = 3.1082368759020804e-16
T covariance error = 1.6051159987815854e-16
```

All-edge sparse and dense unequal-width block matrices are identical to numerical
precision in the canonical validation:

```text
S error = 0.0
H error = 0.0
T error = 0.0
```

## Self-consistent representation-neutral nuclear dynamics

Time-dependent complex-gauge equivalence:

```text
dt=0.0100  coefficient error=1.946455377220016e-12
dt=0.0050  coefficient error=4.858956162973289e-13
dt=0.0025  coefficient error=1.2184279933242558e-13
```

Observed orders:

```text
2.002130949292044
1.995625423493078
```

At the finest step:

```text
position error = 0.0
momentum error = 1.1443916996305594e-16
base generalized norm drift = 4.085620730620576e-14
gauge generalized norm drift = 1.0880185641326534e-14
```

This validates coefficient-coupled nuclear guidance under a genuinely complex,
coordinate-dependent electronic gauge.

## Adaptive block lifecycle

```text
zero-block birth represented-state change = 0.0
zero-block prune projection loss = 0.0
pruned coefficient error = 0.0
```

The insertion rule is exact at the birth event and the pruning reference uses the
nonorthogonal metric Schur complement.

## Generic electronic observable

```text
base expectation = 0.010686571904291245
gauge expectation = 0.010686571904291247
absolute error = 1.734723475976807e-18
```

The validation deliberately contains nonzero imaginary electronic data:

```text
max |Im H| = 0.016864365121461722
max |Im dH| = 0.006399344322890202
max |Im D| = 0.09364391771472179
```

## Full-subspace provider diagnostics

```text
subspace checks = 8
subspace ambiguities = 0
minimum singular value = 0.9999999999999991
```

The provider records Procrustes alignment diagnostics without forcing a local gauge
rotation, so no uncomputed gauge-derivative term is silently inserted into the
connection.

## Complex dtype audit

```text
core files scanned = 9
unclassified suspicious real casts = 0
passed = True
```

Three real casts are intentionally retained only in the inherited spin-free adiabatic
source adapter for energies, gradients, and real NAC input. They are converted into the
complex v0.21 operator contract immediately afterward.

## PySCF status

```text
PySCF installed in build environment: False
real PySCF v0.21.2 runtime trajectory validated: False
SOC Hamiltonian introduced: False
```

A real molecular backend remains an empirical validation milestone before an ab-initio
SOC claim. It is not required before the first exactly reproducible analytic SOC model.

## Claim boundary

v0.21.2 is a pre-SOC-hardened research framework. It does not claim production AIMS,
physical SOC dynamics, or real PySCF runtime validation.
