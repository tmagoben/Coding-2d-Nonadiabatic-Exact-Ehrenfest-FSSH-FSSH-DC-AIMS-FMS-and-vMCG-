# v0.24.1 build validation

Release date: 2026-08-24

## Canonical environment

- CPython 3.12.13 on Linux x86-64, little endian
- NumPy 2.5.2
- SciPy 1.18.0
- h5py 3.16.0
- PySCF 2.13.1 from the preserved hash-locked wheel runtime
- BLAS/OpenMP-related thread controls fixed to one
- OpenMolcas: not installed and not executed

## Recorded validation

- Focused v0.24.1 SOC/release/packaging set: **15 passed in 46.22 s**.
- Canonical acceptance campaign: **315/315 gates passed** (256 inherited, 59 new).
- Real PySCF static-SOC runtime audit: **39/39 gates passed**.
- Spin-algebra and fail-closed controls: **20/20 gates passed**.
- Full cumulative source-tree suite: **418 passed in 429.72 s**.
- Two independent source-safe campaign writes were byte-identical.
- Canonical runtime-evidence SHA-256:
  `9ee704ef66c98b8057630ac6f345dfd42596d6ae20166ad0752b92b740cb2c0f`.
- Canonical campaign SHA-256:
  `91d546b403734b09cbaf728db8832cf39e4a84c46801279e46f899e0dde8f6bd`.
- Two fixed-epoch offline wheel builds were byte-identical; final wheel SHA-256:
  `88eb7a9482ada00a6990c52a991c7809e5c1ffead402b4422a3aac3e4bfb3b01`.
- The final wheel was extracted outside the source tree, imported as v0.24.1, and
  reran the real PySCF runtime audit at **39/39** gates with
  `||H_soc||_F = 139.11503628860513 cm^-1`.
- The clean source candidate contains **858 files** before archiving, with no build,
  egg-info, pytest-cache, bytecode, or `__pycache__` artifacts.

## Claim boundary

This build validates a direct fixed-geometry PySCF BP-SOMF state-interaction SOC
implementation and real OH/STO-3G execution. It does not validate physical SOC
derivatives, derivative connections, cross-geometry many-electron overlaps, trajectory
state tracking, basis/method convergence, general ab-initio accuracy, or a live
trajectory-ready molecular-SOC backend. OpenMolcas external admission remains closed.
