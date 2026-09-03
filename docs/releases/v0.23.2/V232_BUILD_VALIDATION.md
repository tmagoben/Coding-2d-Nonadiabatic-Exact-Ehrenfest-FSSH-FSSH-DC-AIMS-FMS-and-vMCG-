# v0.23.2 build validation

Build date: 2026-08-21

## Source and metadata

- Package version, project metadata, and citation version: 0.23.2.
- Optional PySCF dependency: exactly 2.13.1.
- Exact CPython 3.12/Linux x86-64 runtime wheels: SHA-256 locked.
- Source tree contains no virtual environment, pytest cache, bytecode cache, or
  package egg-info directory at packaging time.
- Python sources pass an AST syntax parse.
- The previous v0.23.1 archive remains sealed and unchanged.

## Test results

The full cumulative test suite was executed with PySCF 2.13.1 and all numerical
thread environments fixed to one:

```text
383 passed in 279.48s (0:04:39)
```

The v0.23.2 packaging, overlap, admission, runtime-contract, and inherited PySCF
adapter focus set separately reported:

```text
29 passed in 1.79s
```

## Canonical release campaign

The real-runtime release campaign passes 168/168 native Boolean gates:

| Gate family | Count | Result |
|---|---:|---|
| Inherited v0.23.1 | 123 | pass |
| Real PySCF runtime | 28 | pass |
| Other v0.23.2 overlap, NAC, and admission controls | 17 | pass |
| **Total** | **168** | **pass** |

The validated positive claims are real PySCF spin-free execution, SA-CASSCF
analytic gradients, NAC/overlap consistency, finite-manifold overlap
contractions, trust-anchored admission logic, and inherited analytic physical-SOC
fixtures. External/live molecular-SOC admission, ab-initio SOC, and live PySCF
SOC remain false.

## Reproducibility

Two consecutive source-tree campaign generations and one installed-only archive
generation produced byte-identical files:

```text
campaign JSON       269e1cd324f5a2050fbc6bfd897051f510d4e95ded59ec4454d1accfc8cce89b
runtime JSON        55e570af5d69b9b2ac63b729542ae76dba5dc8466965c49670160b441164be16
evidence payload    281474f0478cfe6bea2e25f4649e0c05d8b4cfbf407af54d47c959a9c10335fd
runtime environment 3dd0015bd88422905481337ad692f68fd7415ec033fbb44e0d23777ca76b94a8
```

The evidence verifies 1,193 RECORD-hashed PySCF files totaling 168,868,242 bytes.

## Installed-only archive check

The candidate ZIP was extracted to a temporary directory, built as a regular
wheel, and force-installed without dependencies into the pinned environment. An
import from outside the extracted tree resolved to `site-packages` at version
0.23.2. Four focused v0.23.2 test modules then reported:

```text
27 passed in 2.76s
```

The installed-only campaign passed 168/168 and reproduced the canonical JSON
files byte-for-byte. The final ZIP is integrity-tested with `unzip -t`; its
SHA-256 is recorded in the adjacent `.zip.sha256` sidecar.
