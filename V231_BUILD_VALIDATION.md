# v0.23.1 build validation

Build date: 2026-08-21

## Baseline integrity

v0.23.1 was staged independently from the sealed v0.23.0 archive:

```text
Gaussian-Nonadiabatic-Dynamics-v0.23.0.zip
SHA-256 fc39785fc68a2274d67663ee621c42982f1312c15667eb35e444382e5b1983e1
```

The checksum and ZIP integrity were reverified before extraction. The v0.23.0 archive
was not modified.

## Environment and metadata

- Python package version: `0.23.1`
- optional PySCF target: exactly `2.13.1`
- local PySCF/Psi4/xTB runtime: unavailable
- live method-specific SOC adapter: unavailable
- external or live backend admitted: false

No electronic-structure runtime, raw output, or physical accuracy result was mocked.
The canonical datasets are explicitly classified as validation fixtures.

## Focused validation

The complete v0.23.1 focused release set reports:

```text
18 passed in 42.63s
```

The canonical campaign reports:

```text
123/123 gates passed
93 inherited v0.23.0 gates
30 new v0.23.1 gates
```

Reference fingerprints:

```text
even replay   2cabad829ff634ca5612c1e9e8acfbec0afe0a3f6c473ff951fe1fcf73daa651
even dossier  f7d1c3479911ac2b864a12701397eeb99100f0e2038eb4b5a2d0b69b4ffbc82d
odd replay    7cef09125fd86aa0970ee965dcf0ed19448d1d4b91c3fd86341f390d09086797
odd dossier   865d7e88d62b937bef9aee8faf0ffee64031717524d431a3ccef9c561b612487
```

## Cumulative and static validation

All package, test, and example Python sources parse successfully:

```text
425 Python files parsed
```

The cumulative source suite, including every inherited numerical regression and the
v0.23.1 packaging assertion, reports:

```text
353 passed in 240.30s (0:04:00)
```

The machine-readable claim scan found no positive assertion for external snapshot
admission, live backend admission, ab-initio SOC validation, or live PySCF validation.

## Archive validation

The source tree was archived without `.venv`, `.pytest_cache`, `__pycache__`, or
egg-info entries. The archive contains 789 entries and passes `unzip -t` without
compressed-data errors.

The archive was unpacked into a fresh temporary directory and installed with package
indexes disabled. The extracted source package was moved outside the test root before
execution, and the import path confirmed use of the installed site-packages copy.

```text
installed module version:        0.23.1
installed distribution version:  0.23.1
isolated v0.23.1 tests:           18 passed in 41.96s
```

The installed-package example recomputed all 123 gates, both replay files, both
dossiers, and every raw receipt artifact. The complete v0.23.1 artifact set had the
same aggregate hash before and after recomputation:

```text
e29b2a7b192aaef3d6aee02d25c4bd16da193cc95e2bcd62f89429119f5a3afb
```

The final archive SHA-256 is stored in the adjacent `.sha256` sidecar to avoid a
self-referential checksum inside the archive.

## Scientific claim

v0.23.1 validates raw-evidence dossier and calculation-receipt machinery. It does not
admit an external snapshot or live molecular SOC backend and does not validate
ab-initio or PySCF SOC accuracy.
