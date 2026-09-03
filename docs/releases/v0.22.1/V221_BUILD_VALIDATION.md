# v0.22.1 build validation

Date: 2026-08-20

## Build identity

```text
project: Gaussian Nonadiabatic Dynamics
release: v0.22.1
theme: corrective SOC derivative, symmetry, exact-grid, and convergence hardening
physical analytic SOC retained: True
molecular SOC backend admitted: False
ab-initio SOC validated: False
real PySCF runtime validated: False
external magnetic field: False
```

v0.22.1 was staged by extracting the validated v0.22.0 archive into a separate source
tree. The v0.22.0 archive was not modified. Its SHA-256 remained
`206e7555e117e4683f53e2c33cabc42a306178ca0850b957276423ae5020e2cc`, and
`unzip -tq` continued to report no compressed-data errors.

## Source environment

The staged source is installed from its own `pyproject.toml` and reports distribution
and module version 0.22.1. The package resolves to this staged source tree. The build
uses Python 3.12.13, NumPy 2.3.5, SciPy 1.17.0, and pytest 9.1.1. PySCF is not installed
and no PySCF runtime claim is made.

## Focused v0.22.x suite

Command:

```bash
python -m pytest -q tests/test_v0220*.py \
  tests/test_v0221_contract_and_grid.py \
  tests/test_v0221_convergence_and_release.py
```

Result:

```text
38 passed in 97.27s
```

This jointly verifies the inherited physical analytic-SOC behavior and every new
behavioral correction.

## Cumulative source suite

Before adding the packaging-record test to the measured build record:

```text
317 passed in 162.59s
```

The final cumulative run, including the metadata/build-record test, reports:

```text
318 passed in 158.94s
```

## Syntax and claim audit

Every Python source outside generated environments, caches, build directories, and
editable metadata was parsed with `ast.parse`:

```text
Python files parsed: 407
parse failures:      0
```

The canonical JSON contains 67 native JSON booleans and all are true. A source and
documentation scan found zero positive claims that ab-initio SOC, a molecular SOC
backend, or a PySCF SOC runtime had been validated.

## Canonical campaign

Command:

```bash
python examples/120_recompute_v0221_campaign.py
```

Result:

```text
67/67 checks passed
inherited v0.22.0 gates: 53/53
new corrective gates: 14/14
analytic models only: True
molecular SOC backend admitted: False
PySCF runtime validated: False
```

Selected measured results:

```text
generic spin-free component derivative residual: 4.0135369935227635e-14
generic SOC component derivative residual:       2.5196419054680800e-15
cancelled spin-free component residual:           6.6868159521187130e-04
cancelled SOC component residual:                 6.6868159521183930e-04
nonunitary time-reversal residual:                8.4016805041680580e-01
Gaussian-basis fine population difference:        8.9321957644692940e-10
Gaussian-basis narrowing ratio:                   4.8276881041709860e-03
finest sparse-threshold coefficient error:        8.6736174051275750e-19
```

The machine-readable record is
`results/v0221_corrective_hardening_campaign.json`.

## Archive verification

The final archive contains one top-level
`Gaussian-Nonadiabatic-Dynamics-v0.22.1` directory and excludes virtual environments,
Python/test caches, editable-install metadata, bytecode, build/distribution outputs,
editor metadata, operating-system metadata, and VCS paths.

```text
archive entries:                 659
top-level roots:                 1
forbidden archive paths:         0
compressed-data errors:          0
isolated wheel version:          0.22.1
isolated v0.22.1 tests:          17 passed in 56.11s
isolated campaign gates:         67/67
isolated campaign JSON match:    byte-for-byte
```

The isolated import resolved from the fresh installed wheel's `site-packages`, not from
the staging source. The final SHA-256 is recorded in the adjacent `.sha256` sidecar.

## Claim boundary

This build validates corrective contracts around physical analytic singlet–triplet and
Kramers-doublet SOC. It does not validate a molecular or ab-initio SOC backend, a real
PySCF SOC trajectory, external magnetic fields, production AIMS equations, or molecular
predictive accuracy.
