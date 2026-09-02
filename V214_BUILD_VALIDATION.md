# v0.21.4 Build Validation

Date: 2026-08-20

## Build identity

```text
project: Gaussian Nonadiabatic Dynamics
release: v0.21.4
theme: differential-provider and deterministic-restart certification
physical SOC Hamiltonian introduced: False
physical SOC derivative introduced: False
real PySCF runtime validated: False
```

v0.21.4 was staged by extracting the validated v0.21.3 archive into a separate source
tree. The v0.21.3 source archive was not modified; its SHA-256 remained
`21bc8314bdfc2ed2a6d2839b4e5ba6907562e86326bbad44619efd39a0e20456`.

## Editable source install

The staged source was installed from its own `pyproject.toml` and reported:

```text
Python:          3.12.13
gaussian-nadyn:  0.21.4
module version:  0.21.4
NumPy:           2.3.5
SciPy:           1.17.0
pytest:          9.1.1
PySCF installed: False
```

The package path resolved to the v0.21.4 source tree. Explicit setuptools package
discovery remained active.

## Cumulative source suite

Command:

```bash
python -m pytest -q
```

Result:

```text
280 passed in 81.22s
```

This includes all inherited tests and 18 focused v0.21.4 tests covering fixed and
moving-frame provider differentials, wrong-K/wrong-D negative controls, provider
fingerprint admission, dense/sparse/moving-frame restart, integrity and identity
failures, byte manifests, canonical UID edges, retained zero-block densities, adaptive
global-step lifecycle, exact zero-SOC equivalence, metadata, and the release campaign.

## Canonical campaign

Command:

```bash
python examples/118_recompute_v0214_campaign.py
```

Result:

```text
21/21 checks passed
inherited v0.21.3 acceptance: True
physical SOC Hamiltonian introduced: False
physical SOC derivative introduced: False
PySCF runtime validated: False
```

Selected results:

```text
fixed-frame K differential error:        8.232012296296499e-15
moving-frame K differential error:       6.566733600790935e-14
moving-frame D differential error:       1.6631553311401476e-11
maximum overlap-isometry residual:       2.4405271014360246e-16
zero-SOC position/momentum error:         0.0 / 0.0
zero-SOC coefficient error:              7.255325009694841e-17
dense restart coefficient error:         8.238824461827495e-16
sparse restart coefficient error:        8.238824461827495e-16
moving-frame restart coefficient error:  1.9932132794842223e-15
```

The machine-readable result is
`results/v0214_differential_restart_campaign.json`.

## Syntax and claim audit

Every Python source outside generated environments, caches, build directories, and
editable metadata was parsed with `ast.parse`:

```text
Python files parsed: 394
parse failures:      0
```

A source/document/result scan found no claim that physical SOC, physical SOC
derivatives, or a real PySCF runtime had been validated.

## Archive policy

The release archive contains one top-level
`Gaussian-Nonadiabatic-Dynamics-v0.21.4` directory. It excludes virtual environments,
Python/test caches, editable-install metadata, bytecode, editor metadata, and operating
system metadata. The final archive is re-extracted into a new directory for an isolated
editable install, a focused v0.21.4 suite, campaign verification, and archive-path
hygiene checks. Its SHA-256 is recorded in the adjacent `.sha256` sidecar.

The staged archive contained 635 entries and exactly one top-level directory. Its
forbidden-path scan returned no environment, cache, editable metadata, bytecode, build,
distribution, editor, operating-system, or VCS paths. `unzip -tq` reported no compressed
data errors.

After extraction, a new environment installed the package from the extracted
`pyproject.toml`; distribution and module versions both reported 0.21.4 and the import
path resolved inside the extracted tree. With the declared pytest dev dependency
installed, the focused archive suite reported:

```text
18 passed
```

The extracted canonical campaign separately reported 21/21 gates passing, including
inherited v0.21.3 acceptance and the explicit no-physical-SOC/no-PySCF-runtime boundary.

## Claim boundary

This build validates analytic provider consistency, zero-SOC integration plumbing, and
deterministic state continuation. It does not validate a physical SOC Hamiltonian,
physical SOC derivatives, exact-grid SOC dynamics, ab-initio SOC, a real PySCF
trajectory, production AIMS equations, or production asynchronous scheduling.
