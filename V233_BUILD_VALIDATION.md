# v0.23.3 build validation

Release date: 2026-08-21

## Canonical environment

- CPython 3.12.13 on Linux x86-64, little endian
- NumPy 2.5.2
- SciPy 1.18.0
- h5py 3.16.0
- PySCF 2.13.1 from the hash-locked wheel runtime
- BLAS/OpenMP-related thread controls fixed to one

## Recorded validation

- Focused regression after the gauge-optimizer/physical-overlap separation:
  23 passed.
- Full cumulative source-tree suite: **403 passed in 309.15 s**.
- Canonical acceptance campaign: 208/208 gates passed (168 inherited, 40 new).
- Two independent campaign writes were byte-identical; canonical campaign
  SHA-256: `cb0eb5814ad47441b620fe394b631fa01996e94313df79da837a0d88d07980cf`.
- Final source archive: public version/import smoke test and the full 208-gate
  campaign passed from an isolated, non-source installation.

The archive checksum sidecar is verified during release sealing. Generated
Python/test caches and package metadata are excluded from the archive.

## Claim boundary

The validated build includes real spin-free PySCF SA-CASSCF evidence and analytic
physical-SOC fixtures. It does not admit an external or live molecular-SOC source,
does not validate ab-initio SOC accuracy, and does not claim a live PySCF SOC
runtime.
