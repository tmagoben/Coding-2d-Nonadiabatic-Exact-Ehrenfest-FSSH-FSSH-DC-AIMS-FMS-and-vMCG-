# v0.25.0 build validation

Release date: 2026-08-24

## Canonical environment

- CPython 3.12.13 on Linux x86-64, little endian
- NumPy 2.5.2
- SciPy 1.18.0
- h5py 3.16.0
- PySCF 2.13.1 from the preserved hash-locked wheel runtime
- BLAS/OpenMP-related thread controls fixed to one
- OpenMolcas: not installed and not executed

## Scientific acceptance

- Canonical acceptance campaign: **460/460 gates passed**.
- Inherited v0.24.2 gates: **400/400 passed**.
- New deterministic variational-SOC audit: **45/45 passed**.
- New adversarial/core controls: **15/15 passed**.
- Even/odd maximum electronic-norm drift: below `1.8e-14`.
- Maximum signed-reversal spinor residual: below `1.3e-14`.
- Complex-gauge spinor covariance residual: below `1.5e-14`.
- Maximum polar-transport unitarity residual: below `1.6e-15`.
- Canonical evidence JSON SHA-256:
  `b8305a63207a807b0cdcb3bc4e4b7c84aa4334bea8d9e30f7a32e7fcb11bad72`.
- Canonical campaign JSON SHA-256:
  `1e1b9c95357a015d041d13a2631f8f07307ce26dd3fc2bf2e4dafa27694c331a`.
- Full cumulative source-tree suite: **454 passed in 546.93 s**.

## Distribution validation

- Two independent campaign/evidence rewrites were byte-identical.
- Two fixed-epoch offline wheel builds through `setuptools.build_meta` were
  byte-identical; final wheel SHA-256:
  `9525baf972ccf2f67d6a0dcdebc34953fe64151bf86d9830000ed74628819cc0`.
- The wheel was extracted outside the source tree and imported as v0.25.0 without
  PySCF; the complete v0.25.0 API was present.
- That isolated wheel reran the deterministic variational-SOC audit at **45/45
  gates** with evidence fingerprint
  `a598a3f9382e544b96059054f10c31df4166dc3fd7f6d9cc491ae8a6b6094b2c`.
- Wheel metadata reports version 0.25.0, Python >=3.10, NumPy/SciPy core
  dependencies, and the exact optional `pyscf==2.13.1` dependency.
- The final clean source candidate contains **896 files**, with no build, dist,
  egg-info, pytest-cache, bytecode, or `__pycache__` artifacts.
- Two sorted fixed-metadata source archives were built independently and verified
  byte-identical; each was extracted and checked against its source manifest.

## Claim boundary

This build validates the restricted single-packet TDVP integrator on analytic even
and odd SOC models. It does not validate full multi-Gaussian/adaptive-width TDVP,
general noncanonical use of Verlet, a trajectory-ready PySCF molecular-SOC backend,
or general ab-initio SOC dynamics accuracy.
