# v0.23.0 release notes

v0.23.0 adds the molecular spin-orbit-coupling (SOC) backend admission layer. It
does not admit a real molecular backend. Instead, it freezes the capability,
identity, evidence, replay, and failure contracts that a live or externally captured
ab-initio provider must satisfy before it can drive moving-nuclear dynamics.

## Added

- Two explicit capability tiers: `static_soc` and `trajectory_ready`.
- A fingerprinted molecular identity containing the backend and method versions,
  charge and electron count, SOC and scalar-relativistic operators, derivative method,
  units, active space, coordinate definition, state-tracking policy, and—for real
  sources—atomic identities, isotope masses, reference geometry, calculation-input
  hash, and environment hash.
- Independent admission evidence for reference agreement, basis convergence, method
  convergence, translational and rotational invariance, and overlap-based state
  tracking.
- A deterministic manifest/NPZ replay format containing component-resolved operators,
  connections, mass matrices, all pair overlaps, convergence flags, time reversal,
  physical projectors, and complete provenance.
- Separate protocol and real-backend outcomes. An analytic fixture may pass the former
  but can never pass the latter.
- A fail-closed PySCF boundary that requires an installed runtime and an injected,
  method-specific provider with affirmative SCF, correlated, and SOC convergence.
- Even-electron singlet–triplet and odd-electron two-Kramers-doublet replay fixtures.
- Five independent evidence-negative controls plus corruption, unconverged-record,
  coordinate-miss, cross-dataset, and static-only controls.

## Acceptance

The canonical campaign passes **93/93 gates**: all 67 v0.22.1 gates and 26 new
v0.23.0 gates. The deterministic fixtures reproduce every stored operator component
and overlap exactly.

Reference replay fingerprints:

- even singlet–triplet:
  `a0c90420ace96c899b5e033a8b03d43cd316eb17511985140541be4e4dec8255`
- odd two-doublet:
  `9e12baf2950ee2f7fe6ec5debe183cc3c5b2ff1e25ef3df20fe636e05a2c9acb`

Recompute the campaign with:

```bash
python examples/121_recompute_v0230_campaign.py
```

## Claim boundary

v0.23.0 validates the molecular-SOC admission protocol and deterministic replay
machinery. It does **not** validate ab-initio SOC accuracy, admit a real molecular SOC
backend, or validate a live PySCF SOC runtime. PySCF was not installed in the release
environment. v0.23.1 is the first possible live-backend admission release, conditional
on a method-specific implementation and complete independent evidence.
