# v0.23.1 PySCF adapter boundary

> **Historical convention notice:** v0.23.2 real-runtime finite-difference
> evidence supersedes the tuple interpretation below. Production internal
> `d[i,j]` now uses PySCF `state=(i,j)`, `mult_ediff=False`, and
> `use_etfs=False`. See `V232_NAC_CONVENTION_ERRATUM.md`.

## Frozen runtime

The optional release dependency is pinned to PySCF 2.13.1. PySCF was not installed in
the build environment, so this pin is an integration target rather than a validated
runtime result.

The official SA-CASSCF NAC API defines `state=(ket,bra)` and returns
$\langle\mathrm{bra}|\partial_R\mathrm{ket}\rangle$. v0.23.1 freezes that orientation
as a literal adapter convention. See the
[official PySCF NAC API](https://pyscf.org/pyscf_api_docs/pyscf.nac.html).

## What core importability does not prove

Importing PySCF does not establish:

- a chosen state-interaction SOC Hamiltonian and prefactors;
- physical nuclear derivatives of that SOC Hamiltonian;
- complete even- or odd-electron multiplet construction;
- many-electron overlaps across geometries;
- raw-output parsing and replay reconstruction;
- molecular reference, basis, method, frame, or tracking accuracy.

The framework deliberately does not infer those features from AO integrals, X2C
availability, or an SCF object.

## Required method-specific engine

`PySCFMethodSpecificSOCAdapterV231` accepts only an injected engine that declares and
implements all of:

1. state-interaction SOC;
2. physical SOC derivatives;
3. analytic spin-free gradients;
4. derivative connections;
5. many-electron cross-geometry overlaps;
6. a raw-artifact parser;
7. a fresh execution in the pinned runtime.

It must expose component operators, snapshots, overlaps, raw-artifact writing, and
raw-artifact validation. Every emitted snapshot must affirm SCF, correlated-state, SOC,
derivative, and overlap convergence.

## Admission sequence

```text
installed pinned runtime
  -> complete method-specific engine
  -> fresh calculations and raw artifacts
  -> exact replay capture
  -> deterministic dossier
  -> executable artifact validation
  -> v0.23.0 physical/symmetry audit
  -> v0.23.1 external or live admission
```

The sequence stops at the first missing condition. In this release it stops at the
first condition because PySCF is unavailable.
