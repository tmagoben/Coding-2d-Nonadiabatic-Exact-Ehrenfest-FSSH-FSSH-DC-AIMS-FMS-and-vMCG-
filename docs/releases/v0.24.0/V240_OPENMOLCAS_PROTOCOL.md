# v0.24.0 OpenMolcas RASSI-SO protocol

## Status

This document freezes an intake protocol. It is not evidence that OpenMolcas was
executed and is not an ab-initio SOC validation claim.

## Electronic-structure target

- Backend: OpenMolcas 26.06.
- Molecule: neutral H2O, 10 electrons, explicit O-16/H-1 isotope masses.
- Geometry and all displacements: bohr, C1 symmetry.
- Orbital basis: ANO-RCC-VDZP.
- Scalar relativity: second-order Douglas-Kroll-Hess (DKH2).
- Active space: CAS(8,6), held identical by definition across all records.
- Spin-free energies: state-specific CASSCF followed by single-state CASPT2.
- RASSI diagonal energies: `EJOB`.
- SOC operator: RASSI-SO using SEWARD AMFI integrals.
- RASSI selection: `NROFJOBIPH=2 1 1;1;1`, `SPINORBIT`, `EJOB`.
- Ordered frame: `S0(M=0)`, `T1(M=-1)`, `T1(M=0)`, `T1(M=+1)`.
- Magnetic field: zero.

These are project choices. The official OpenMolcas RASSI documentation motivates
the module/keyword flow, but it does not validate this project's selected active
space, basis, molecule, convergence, or accuracy.

Primary implementation references:

- [OpenMolcas RASSI user guide](https://molcas.gitlab.io/OpenMolcas/sphinx/users.guide/programs/rassi.html)
- [OpenMolcas RASSI tutorial](https://molcas.gitlab.io/OpenMolcas/sphinx/tutorials/tut_rassi.html)
- [OpenMolcas SEWARD relativistic operators](https://molcas.gitlab.io/OpenMolcas/sphinx/users.guide/programs/seward.html)

## Displaced-geometry inventory

The bundle has one reference calculation plus `+h` and `-h` calculations for each
of 9 Cartesian coordinates at 0.004, 0.002, and 0.001 bohr: exactly 55 records.
Every record carries separately converged Gateway, Seward, SCF, singlet RASSCF,
triplet RASSCF, singlet CASPT2, triplet CASPT2, and RASSI-SO stages.

Cross-geometry overlaps are produced by an independently identified biorthogonal
CASSCF wavefunction-overlap exporter. They are not inferred from RASSI's
within-geometry overlap matrix. Raw finite-manifold contractions are retained;
operator transport uses their certified right-to-reference unitary polar factors.

## Raw artifact set

Each record directory contains exactly:

```text
openmolcas.input
openmolcas.output
rassi.h5
gnd_rassi_export_v240.json
```

The root contains the strict manifest, an independent validation artifact, and a
content-addressed raw blob for every reference, basis, method, frame, and tracking
observation. SHA-256 digests bind every native artifact to its export, every validation
observation to a raw blob, every export to the exact protocol, and the root manifest
to an out-of-band admission policy.

The checked protocol fixture deliberately uses a non-HDF5 placeholder and a mandatory
`NO-OPENMOLCAS-EXECUTION` marker. Relabeling it as external evidence is rejected.
