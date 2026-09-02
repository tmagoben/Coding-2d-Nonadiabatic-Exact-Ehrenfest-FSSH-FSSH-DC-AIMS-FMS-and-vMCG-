# v0.23.2 release notes

v0.23.2 validates the pinned PySCF 2.13.1 spin-free runtime, corrects the
SA-CASSCF nonadiabatic-coupling (NAC) adapter orientation using many-electron
overlap finite differences, replaces an invalid cross-geometry isometry demand
with the finite-manifold contraction contract, and hardens the molecular-SOC
admission boundary. It does **not** admit or validate an ab-initio SOC backend.

## Corrected

- The production PySCF adapter now requests `state=(i,j)` for the internal
  `d[i,j]=<Phi_i|d Phi_j/dR>` convention. In PySCF 2.13.1 this is the mapping
  observed from phase-aligned many-electron overlap derivatives. The earlier
  literal interpretation `state=(j,i)` had the opposite sign.
- Full overlap derivatives use `use_etfs=False`. ETF-corrected NACs remove the
  translational component and are retained as a distinct diagnostic quantity.
- Cross-geometry overlaps in a selected finite state space need not be unitary.
  Self overlaps must be identity, reverse overlaps must be adjoints, and every
  cross-overlap singular value must be at most one within tolerance. Singular
  values below one represent amplitude leaking into omitted roots.

## Added

- A pinned, hash-locked CPython 3.12/Linux x86-64 PySCF runtime specification.
- Exact module and distribution version checks plus RECORD-content verification
  and reproducible runtime fingerprints.
- Real H3+ SA-CASSCF(2e,3o)/STO-3G evidence for three equally weighted singlet
  roots: energies, analytic gradients, NACs, many-electron overlaps, and three
  central-difference steps.
- Overlap diagnostics reporting self-identity, reciprocity, singular-value
  retention, contraction excess, and isometry defect separately.
- Structure-before-capability checks, callable-interface enforcement, exact
  method identity, seven convergence stages, typed parser proof, trust-anchored
  validator identity, and unmergeable live/external admission evidence.

## Acceptance and claim boundary

The canonical campaign passes **168/168 gates**: 123 inherited gates, 28 real
PySCF runtime gates, and 17 other v0.23.2 correction/admission gates (45 new
gates total). It validates the spin-free PySCF runtime, SA-CASSCF gradients,
NAC/overlap consistency, the finite-manifold overlap contract, and the
trust-anchored admission machinery.

External molecular-SOC snapshots, live molecular-SOC backends, ab-initio SOC,
and a live PySCF SOC runtime remain explicitly **not admitted**.

Recompute with:

```bash
python examples/123_recompute_v0232_campaign.py
```
