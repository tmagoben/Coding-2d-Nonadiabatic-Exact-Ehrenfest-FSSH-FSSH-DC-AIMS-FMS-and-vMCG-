# v0.24.0 release notes

v0.24.0 introduces the first method-specific external molecular-SOC intake path. It
freezes an OpenMolcas RASSI-SO protocol, validates all parsing/admission machinery with
a conspicuously synthetic protocol fixture, and keeps the external claim closed until
real artifacts are supplied.

## Added

- Exact OpenMolcas 26.06 H2O CAS(8,6)/CASPT2/RASSI-SO/AMFI protocol identity.
- One reference plus 54 centered Cartesian displacement records.
- Strict input, output, HDF5, export, validation-blob, and manifest digest binding.
- Explicit fixture versus external source classes and synthetic-relabel rejection.
- Separate transported spin-free and SOC derivative convergence audits.
- Independent reference, basis, method, rigid-frame, and tracking evidence schema.
- Caller-owned parser, convention, manifest, environment, and exporter trust anchors.
- Admission-bound frozen-snapshot unitary propagation, zero-SOC mode, and restart.
- 48 new gates, for 256 cumulative release gates.

## Claim boundary

The protocol fixture is validated. No OpenMolcas runtime exists in the build
environment, no genuine OpenMolcas output was supplied, and therefore no external or
live molecular-SOC source, ab-initio SOC accuracy, or OpenMolcas execution is claimed.
The independent native HDF5/text numerical cross-parser is not yet implemented, and
the corresponding admission prerequisite is hard-coded false.

Recompute with:

```bash
python examples/126_recompute_v0240_campaign.py
```
