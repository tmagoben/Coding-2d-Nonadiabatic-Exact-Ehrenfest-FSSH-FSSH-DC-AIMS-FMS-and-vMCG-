# v0.24.2 build validation

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

- Focused v0.24.2 core/runtime/packaging set: **15 passed in 1.71 s**.
- Canonical acceptance campaign: **400/400 gates passed**.
- Inherited v0.24.1 gates: **315/315 passed**.
- Real connected-geometry PySCF audit: **60/60 passed**.
- New adversarial/core controls: **25/25 passed**.
- Direct-JK/explicit maximum absolute difference:
  `1.7763568394002505e-15`.
- Minimum retained endpoint-overlap singular value: `0.9968178663528345`.
- Maximum matrix/symmetry residual: `5.970344688601405e-12`.
- Canonical evidence JSON SHA-256:
  `43f62aa3b8a04ff7a4467f44ab21d5c914a4130760552633ac8a892206ef8f41`.
- Canonical campaign JSON SHA-256:
  `5b1af33354141b5b432733cda56dfb8a0e6dbdfc252abd444dc42ef2723e19e2`.
- Full cumulative source-tree suite: **434 passed in 490.00 s**.

## Distribution validation

- Two independent real campaign/evidence rewrites were byte-identical.
- Two fixed-epoch offline wheel builds through `setuptools.build_meta` were
  byte-identical; final wheel SHA-256:
  `c57c4a8578fbe5314d8e07ab463503efcb6e415ff9cb78f77454f8b07a7f1306`.
- The wheel was extracted outside the source tree and imported as v0.24.2 without
  PySCF; the optional runtime probe correctly remained unavailable while the public
  v0.24.2 API was present.
- The same extracted wheel, with only the pinned PySCF site-packages added, reran the
  real connected-geometry audit at **60/60 gates**, retained six endpoint receipts,
  and reproduced the `1.7763568394002505e-15` direct-JK/oracle error.
- Wheel metadata reports version 0.24.2, Python >=3.10, NumPy/SciPy core
  dependencies, and the exact optional `pyscf==2.13.1` dependency.
- The final clean source candidate contains **876 files**, with no build, dist,
  egg-info, pytest-cache, bytecode, or `__pycache__` artifacts.
- Two sorted fixed-metadata source archives were built independently and verified
  byte-identical; each was extracted and checked against its source manifest.

## Claim boundary

This build validates direct-JK SOMF execution, connected OH snapshot/overlap evidence,
complete-doublet polar transport, and one-coordinate transported finite differences.
It does not validate continuous physical SOC-state derivative connections, full
Cartesian or analytic derivatives, a real mixed-multiplicity case, general ab-initio
accuracy, or a trajectory-ready/live molecular-SOC backend.
