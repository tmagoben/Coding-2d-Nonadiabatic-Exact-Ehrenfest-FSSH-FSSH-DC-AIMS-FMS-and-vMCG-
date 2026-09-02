# v0.22.0 build validation

Date: 2026-08-20

## Build identity

```text
project: Gaussian Nonadiabatic Dynamics
release: v0.22.0
theme: first physical analytic SOC with singlet–triplet and Kramers-doublet references
physical SOC Hamiltonian introduced: True
physical SOC derivative introduced: True
analytic models only: True
ab-initio SOC validated: False
real PySCF runtime validated: False
external magnetic field: False
```

v0.22.0 was staged by extracting the validated v0.21.4 archive into a separate source
tree. The v0.21.4 archive was not modified. Its SHA-256 remained
`1c511c68fc8a3568dbb21268d5ba402ed84aa026ffb5e6da94b4149166174501`, and
`unzip -tq` continued to report no compressed-data errors.

## Editable source install

The staged source was installed from its own `pyproject.toml` and reported:

```text
Python:          3.12.13
gaussian-nadyn:  0.22.0
module version:  0.22.0
NumPy:           2.3.5
SciPy:           1.17.0
pytest:          9.1.1
PySCF installed: False
```

Explicit setuptools package discovery remained active, and the package resolved to the
staged v0.22.0 source tree.

## Focused release suite

Command:

```bash
python -m pytest -q tests/test_v0220*.py
```

Result:

```text
23 passed in 37.10s
```

The focused suite covers operator composition, analytic physical derivatives, time
reversal, Kramers pairing, gauge-transformed antiunitary representations and projectors,
wrong-derivative and broken-symmetry controls, zero-SOC equivalence, exact-grid
propagation and convergence, Gaussian/grid population agreement, SOC-active restart,
release metadata, and the canonical campaign.

## Cumulative source suite

Command:

```bash
python -m pytest -q
```

Result:

```text
302 passed in 124.68s
```

This contains the full inherited regression history plus the v0.22.0 tests.

## Canonical campaign

Command:

```bash
python examples/119_recompute_v0220_campaign.py
```

Result:

```text
53/53 checks passed
inherited v0.21.4 gates: 21/21
new physical-SOC gates: 32/32
analytic models only: True
ab-initio SOC validated: False
PySCF runtime validated: False
```

Selected numerical results:

```text
singlet–triplet K differential error:       1.8665508646787508e-15
singlet–triplet SOC-force error:            3.108180623686077e-15
doublet K differential error:               2.8189256543263856e-15
maximum Kramers pair splitting:             4.336808689942018e-18
exact-grid observed timestep order:         2.0000059562889683
grid-spacing/box population error:          8.632959798415829e-18
singlet–triplet Gaussian/grid pop. error:   7.405649891539424e-09
doublet Gaussian/grid population error:     2.2499021452878545e-09
dense restart coefficient error:            1.3877787809430044e-17
sparse restart coefficient error:           1.6531559496939515e-16
moving-gauge dynamics coefficient error:    5.936673241766072e-16
```

The machine-readable result is
`results/v0220_physical_analytic_soc_campaign.json`.

## Syntax and claim audit

Every Python source outside generated environments, caches, build directories, and
editable metadata was parsed with `ast.parse`:

```text
Python files parsed: 402
parse failures:      0
```

A source/document scan confirmed that every molecular, ab-initio, PySCF SOC, external
field, and production-AIMS mention remains either an explicit exclusion or a future
milestone. The release claims physical analytic SOC only.

## Archive and isolated-install verification

The final archive contains one top-level
`Gaussian-Nonadiabatic-Dynamics-v0.22.0` directory. It excludes virtual environments,
Python/test caches, editable-install metadata, bytecode, build/distribution outputs,
editor metadata, operating-system metadata, and VCS paths.

```text
archive entries:                 647
forbidden archive paths:         0
compressed-data errors:          0
isolated focused tests:          23 passed
isolated campaign gates:         53/53
isolated campaign JSON match:    byte-for-byte
```

The extracted environment installs distribution and module version 0.22.0, resolves
the import path inside the extracted tree, and independently reproduces the campaign.
The final SHA-256 is recorded in the adjacent `.sha256` sidecar.

## Claim boundary

This build validates physical analytic singlet–triplet and Kramers-doublet SOC,
analytic physical SOC derivatives, exact-grid reference dynamics, complex-gauge
covariance, and deterministic SOC-active restart. It does not validate an ab-initio
or molecular SOC backend, a real PySCF SOC trajectory, external magnetic fields,
production AIMS equations, or molecular predictive accuracy.
