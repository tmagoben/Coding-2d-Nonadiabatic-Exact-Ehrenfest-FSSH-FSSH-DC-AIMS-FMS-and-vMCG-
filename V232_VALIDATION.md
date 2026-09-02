# v0.23.2 validation

## Canonical campaign

`examples/123_recompute_v0232_campaign.py` writes the campaign and standalone
runtime evidence under `results/`. The campaign passes **168/168 gates**:

- 123 inherited v0.23.1 gates;
- 28 real PySCF runtime and NAC/overlap gates;
- 5 finite-manifold overlap controls;
- 10 runtime-admission hardening controls;
- 2 explicit NAC-mapping controls.

The 45 new gates are all native Booleans. The runtime evidence contains exactly
28 passing checks and a canonical SHA-256 digest.

## Real-runtime numerical evidence

- PySCF module and distribution: exactly 2.13.1.
- Real H3+ SA-CASSCF(2e,3o)/STO-3G, three singlet roots.
- Analytic gradients: finite, nontrivial, translation residual below 1e-15.
- NACs: finite, nontrivial, antisymmetry residual at numerical zero.
- Self overlap: identity residual below 2e-15.
- Cross overlap: reciprocity residual below 1e-15 and no contraction excess.
- Production NAC/overlap maximum errors: approximately 7.48e-5, 7.48e-7,
  and 7.48e-9 for tenfold step refinement.
- Scaled NAC tuple symmetry and energy-gap relation: numerical precision.
- ETF and full-overlap derivatives: explicitly demonstrated to be distinct.

## Validated claims

- real PySCF spin-free runtime;
- real PySCF SA-CASSCF analytic gradients;
- real PySCF NAC/many-electron-overlap consistency;
- finite-manifold overlap contraction contract;
- trust-anchored runtime-admission machinery;
- inherited analytic physical-SOC fixtures.

## Explicitly unvalidated

- external molecular-SOC snapshot admission;
- live molecular-SOC backend admission;
- ab-initio SOC accuracy;
- live PySCF SOC runtime.

The final build record is in `V232_BUILD_VALIDATION.md`.
