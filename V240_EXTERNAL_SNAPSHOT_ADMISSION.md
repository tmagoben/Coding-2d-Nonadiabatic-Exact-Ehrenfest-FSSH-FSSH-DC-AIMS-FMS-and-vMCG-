# v0.24.0 external snapshot admission

## Admission sequence

```text
caller-owned trust policy
  |-- exact protocol + SOC convention
  |-- exact manifest/environment digests
  |-- exact parser/exporter identities
  v
strict artifact inventory and SHA-256 binding
  v
native OpenMolcas completion + HDF5 identity
  v
component Hermiticity + time reversal + complete S/T manifold
  v
55-record transported finite-difference derivative evidence
  v
independent reference, basis, method, frame, and tracking evidence
  v
external snapshot admission
  v
admission-bound frozen-nuclei electronic dynamics
```

No bundle can supply its own trust anchor. Parser subclasses are rejected, unknown
files are rejected, paths cannot escape the bundle, and symlinks cannot substitute
for regular artifacts. A source must be explicitly classified as
`external_ab_initio_snapshot`, contain native HDF5 signatures, and carry required
OpenMolcas completion markers.

Protocol validity and physical admission are separate outcomes. The deterministic
fixture passes the former and fails the latter. Consequently v0.24.0 currently has:

- protocol fixture: validated;
- external molecular-SOC snapshot: not admitted;
- live molecular-SOC backend: not admitted;
- ab-initio SOC accuracy: not validated.

The v0.24.0 parser verifies the exact bundle inventory, native-file signatures,
completion markers, and all native/export digests. It does not yet independently
reconstruct the exported SOC matrices from OpenMolcas-native HDF5/text datasets.
`native_numeric_crosscheck` is therefore hard-coded false and is an external-admission
prerequisite. This prevents even a structurally perfect submitted bundle from opening
the gate prematurely.

## Evidence still required to open the gate

1. The complete 55-record OpenMolcas 26.06 archive.
2. A content-addressed runtime/environment record.
3. Native `rassi.h5` and text outputs for every record.
4. The independently identified cross-geometry overlap export.
5. A reference calculation from a backend other than OpenMolcas.
6. Three-level basis and method ladders with raw artifact digests.
7. Translation/rotation calculations and connected-manifold tracking evidence.
8. A reviewer-supplied manifest digest and environment digest for the policy.
9. A version-pinned native OpenMolcas HDF5/text numerical cross-parser that reproduces
   the exported matrices independently.

Only after those artifacts pass the same campaign may the external claim change to
true. A live-backend claim remains a later and independent gate.
