# v0.21.3 Build Validation

Date: 2026-08-13

## Build identity

```text
project: Gaussian Nonadiabatic Dynamics
release: v0.21.3
theme: SOC-contract freeze and degeneracy-safe pre-integration procedures
physical SOC Hamiltonian introduced: False
real PySCF runtime validated: False
```

v0.21.3 was staged as a separate source tree from the validated v0.21.2 archive. The
v0.21.2 source and release ZIP were not modified.

## Isolated editable install

The release was installed into a new virtual environment from its own `pyproject.toml`:

```text
Python:          3.12.13
gaussian-nadyn:  0.21.3
module version:  0.21.3
NumPy:           2.3.5
SciPy:           1.17.0
pytest:          9.1.1
```

The package path resolved to the v0.21.3 source tree and explicit setuptools package
discovery remained active.

## Cumulative test suite

Command:

```bash
python -m pytest -q
```

Result:

```text
263 passed in 47.04s
```

This includes all inherited tests plus 17 v0.21.3 tests covering strict invariants,
model-space/provenance contracts, degeneracy-safe guidance, transactional corrector
state, arbitrary-state projection, complex cache behavior, release metadata, and the
canonical acceptance campaign.

## Canonical campaign

Command:

```bash
python examples/117_recompute_v0213_campaign.py
```

Result:

```text
20/20 checks passed
inherited v0.21.2 acceptance: True
physical SOC Hamiltonian introduced: False
PySCF runtime validated: False
```

Selected numerical results:

```text
rejected Hermiticity defect:       8.944272e-07
current degenerate gauge error:    2.220446049250313e-16
retained degenerate gauge error:   1.1102230246251565e-16
integrated runner norm drift:      5.551115123125783e-16
corrector guide-state rollbacks:   9
four-state projection fidelity:    1.0
projection relative residual:      1.3638610398823859e-27
complex cache round-trip error:    0.0
cached imaginary signal:           0.011539731181531185
```

The machine-readable result is
`results/v0213_soc_contract_freeze_campaign.json`.

## Syntax and source audit

Every Python source outside the virtual environment and generated cache directories was
parsed with `ast.parse`:

```text
Python files parsed: 385
parse failures:      0
```

The build also verifies that real coordinates and mass matrices cannot silently discard
an imaginary component, non-finite provenance cannot be fingerprinted, moving-frame
wavefunction snapshots cannot enter the fixed-frame cache, and a cross-geometry cache
overlap must match the declared fixed frame.

## Archive policy

The release archive contains one top-level
`Gaussian-Nonadiabatic-Dynamics-v0.21.3` directory. It excludes:

- `.venv`;
- `__pycache__` and `.pytest_cache`;
- editable-install `*.egg-info`;
- `.pyc` and `.pyo` files;
- editor and operating-system metadata.

The final ZIP is accompanied by a SHA-256 sidecar and is re-extracted into a separate
directory for focused v0.21.3 tests and campaign verification.

The staged archive contained 619 entries and no forbidden environment, cache,
editable-install metadata, or bytecode paths. Its extracted editable install reported
both package and module version 0.21.3 and resolved imports from the extracted source.
The focused archive suite reported:

```text
17 passed in 7.91s
```

That focused suite includes the 20-gate campaign, so inherited v0.21.2 acceptance and
the explicit no-physical-SOC boundary are both re-evaluated after extraction.

## Claim boundary

This build validates procedures and interfaces for later SOC integration. It does not
validate a physical SOC Hamiltonian, ab-initio SOC derivatives, a real PySCF trajectory,
production AIMS equations, or production asynchronous electronic-structure scheduling.
