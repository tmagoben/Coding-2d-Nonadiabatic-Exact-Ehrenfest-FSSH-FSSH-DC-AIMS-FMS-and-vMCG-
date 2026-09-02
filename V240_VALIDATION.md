# v0.24.0 validation

## Cumulative acceptance

The canonical campaign passes 256/256 native-Boolean gates:

| Layer | Gates |
|---|---:|
| Inherited v0.23.3 | 208 |
| Frozen OpenMolcas protocol | 12 |
| Bundle parser, derivatives, and independent evidence | 16 |
| Fail-closed admission and corruption controls | 12 |
| Frozen-snapshot propagation and restart | 8 |
| **Total** | **256** |

## New numerical checks

- 55 exact record identities: reference plus 9 coordinates, 3 steps, and 2 signs.
- Separate transported `H_spin_free` and `H_SOC` centered differences.
- Hermiticity, time reversal, spin-component degeneracy, finite-manifold retention,
  complete-manifold assignment, and finest-step convergence.
- Independent SOC reference agreement, three-level basis and method ladders, rigid
  translation/rotation spectra, and tracking quality.
- Static matrix-exponential reference, norm preservation, exact zero-SOC path, and
  deterministic checkpoint continuation.

## New negative controls

- fixture-to-external relabel;
- native output corruption;
- digest-consistent native-input/export geometry disagreement;
- unknown files;
- environment and SOC-convention trust-anchor mismatch;
- parser subclass substitution;
- validation JSON and content-addressed raw-blob corruption;
- production dynamics without external admission;
- checkpoint/bundle identity mismatch.

## Deliberately false claims

The campaign asserts that external snapshot admission, live backend admission,
ab-initio SOC validation, OpenMolcas execution, and the native OpenMolcas numerical
cross-parser are all false. These are release invariants, not missing test results.
