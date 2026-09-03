# v0.23.2 pinned PySCF runtime

## Reproducible environment

The validated target is CPython 3.12 on Linux x86-64 with PySCF 2.13.1,
NumPy 2.5.2, SciPy 1.18.0, and h5py 3.16.0. The exact wheel hashes are frozen in
`requirements-pyscf-v232-linux-x86_64-py312.txt`.

```bash
python3.12 -m venv .venv-pyscf-v232
.venv-pyscf-v232/bin/python -m pip install --require-hashes \
  -r requirements-pyscf-v232-linux-x86_64-py312.txt
.venv-pyscf-v232/bin/python -m pip install --no-deps -e .
```

The runtime probe requires both the imported module and installed distribution
to report exactly 2.13.1 and requires the SA-CASSCF NAC API. It fingerprints the
Python executable, PySCF module, NAC module, distribution inventory, every
RECORD-hashed file, scientific dependencies, platform, byte order, and numerical
thread environment. Importability alone is not acceptance.

## Validated calculation

The real-runtime fixture is an asymmetric H3+ triangle in Bohr with charge +1,
spin 0, STO-3G, SA-CASSCF(2e,3o), and three equally weighted singlet roots. Two
nearby geometries establish nonzero energy response, nontrivial analytic
gradients, translational gradient invariance, antisymmetric NACs, identity self
overlaps, reciprocal cross overlaps, and physically contractive selected-state
overlaps.

For the H[2] y coordinate and state pair (0,2), the production derivative is
approximately -0.0312907780 bohr^-1. Central differences at 1e-2, 1e-3, and
1e-4 bohr have maximum errors approximately 7.48e-5, 7.48e-7, and 7.48e-9,
showing the expected second-order convergence.

## Sandbox memory telemetry

Some PID namespaces expose `/proc/self/statm` but not `/proc/<os.getpid()>/statm`.
The guarded release path may substitute the equivalent self-process lookup for
PySCF memory telemetry only. It does not change integrals, tolerances, states,
gradients, NACs, overlaps, or convergence. The canonical evidence explicitly
records `proc_self_statm_requested`.

## Scope

This runtime proves real **spin-free** PySCF SA-CASSCF execution. It does not
provide the method-specific state-interaction SOC Hamiltonian, physical SOC
derivatives, complete spin manifolds, or raw SOC parser needed for molecular-SOC
admission.
