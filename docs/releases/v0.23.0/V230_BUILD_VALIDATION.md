# v0.23.0 build validation

Build date: 2026-08-21

## Baseline integrity

v0.23.0 was staged independently from the immutable v0.22.1 archive:

```text
Gaussian-Nonadiabatic-Dynamics-v0.22.1.zip
SHA-256 25fb0669183bae7dd9628d35fe6802341677a32b6870ab6bf0cff4ccad91cd9c
```

The baseline checksum and ZIP integrity were reverified before staging. The v0.22.1
archive was not modified.

## Environment and metadata

- Python package version: `0.23.0`
- installed distribution version: `0.23.0`
- internal units: bohr, hartree, hartree/bohr
- PySCF import probe: unavailable
- live PySCF SOC adapter validated: false

The unavailable PySCF runtime is an explicit validation result. No dependency was
mocked and no live calculation is claimed.

## Focused validation

The contract, replay/admission, and release campaign tests report:

```text
18 passed in 56.22s
```

The canonical release campaign reports:

```text
93/93 gates passed
67 inherited v0.22.1 gates
26 new v0.23.0 gates
```

Reference replay fingerprints:

```text
even singlet–triplet  a0c90420ace96c899b5e033a8b03d43cd316eb17511985140541be4e4dec8255
odd two-doublet       9e12baf2950ee2f7fe6ec5debe183cc3c5b2ff1e25ef3df20fe636e05a2c9acb
```

Every component, total operator, and cross-record overlap round-trips with exactly zero
reported error. Repeated captures have identical manifest bytes, NPZ bytes, and dataset
fingerprints.

## Cumulative validation

The cumulative source suite, including every inherited numerical regression and the
v0.23.0 packaging assertion, reports:

```text
336 passed in 238.44s (0:03:58)
```

## Static validation

All package, test, and example Python sources parse successfully:

```text
416 Python files parsed
```

A release claim scan must retain false values for real molecular SOC admission,
ab-initio SOC validation, and live PySCF SOC validation.

## Archive validation

The source tree was archived without `.venv`, `.pytest_cache`, `__pycache__`, or
egg-info entries. The archive contains 683 entries and passes `unzip -t` with no
compressed-data errors.

The archive was unpacked into a fresh temporary directory and installed with package
indexes disabled. The extracted source package was then moved out of the test root;
the import path confirmed that tests used the newly installed site-packages copy.

```text
installed module version:       0.23.0
installed distribution version: 0.23.0
isolated v0.23.0 tests:          19 passed in 58.47s
```

The installed-package example recomputed all 93 campaign gates. The campaign JSON and
both replay datasets had the same combined hash before and after recomputation:

```text
208009bce50bddffaa5f807a3fa565d7519f699c339bb183345902c533cd3dd8
```

The final archive SHA-256 is recorded in the adjacent `.sha256` sidecar to avoid a
self-referential checksum inside the archive.

## Scientific claim

v0.23.0 validates the molecular-SOC protocol and deterministic replay implementation.
It admits no real molecular SOC backend and validates neither ab-initio SOC accuracy nor
a live PySCF SOC runtime.
