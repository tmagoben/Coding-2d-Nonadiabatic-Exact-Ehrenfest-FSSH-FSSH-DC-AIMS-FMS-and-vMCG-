# v0.23.3 release notes

v0.23.3 is a corrective compatibility release between the validated spin-free
PySCF foundation and a future method-specific molecular-SOC backend. It makes
finite-manifold transport, replay migration, derivative-coupling identity,
complete multiplet tracking, SOC matrix conventions, and runtime compatibility
explicit and fail closed.

## Added

- A physical finite-manifold overlap object that preserves the raw contraction
  and derives a separate right-to-left unitary polar transport.
- Independent physical-consistency and trajectory-readiness criteria, including
  singular-value retention, condition number, and principal-angle gates.
- Replay format 2 with raw overlaps, unitary transports, singular values, the
  exact overlap policy, and a fingerprinted NAC convention.
- Explicit legacy migration attestations. Unknown or known-wrong NAC signs are
  quarantined; v0.23.3 never silently flips stored data.
- Convention-complete provider identities for replay, cache, and checkpoint use.
- Complete singlet/triplet and Kramers-doublet manifold tracking under independent
  endpoint gauges, including leakage and time-reversal checks.
- A frozen molecular-SOC matrix convention covering operator family, prefactor,
  scalar-relativistic treatment, state order, units, fixed-frame derivatives,
  complete multiplets, and field policy.
- Separate `release_locked` and `scientifically_compatible` runtime profiles.

## Corrected

- State transport no longer treats a nonunitary finite-manifold contraction as
  an electronic coefficient map.
- Procrustes consumers reject rank-lost overlaps instead of returning an
  arbitrary singular-vector completion.
- Legacy replay/cache data without an exact NAC identity no longer enter a
  v0.23.3 trajectory.
- Runtime portability no longer implies byte identity with the canonical build.

## Acceptance and claim boundary

The canonical campaign passes **208/208 gates**: all 168 v0.23.2 gates plus 40
new v0.23.3 controls. The release inherits real PySCF 2.13.1 spin-free
SA-CASSCF validation and physical analytic SOC fixtures.

External molecular-SOC snapshots, live molecular-SOC backends, ab-initio SOC,
and live PySCF SOC remain explicitly **not admitted**.

Recompute with:

```bash
python examples/124_recompute_v0233_campaign.py
```
