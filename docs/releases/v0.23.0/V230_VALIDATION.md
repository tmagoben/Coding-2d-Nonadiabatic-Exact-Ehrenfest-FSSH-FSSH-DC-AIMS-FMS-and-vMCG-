# v0.23.0 validation

## Layered strategy

v0.23.0 keeps implementation checks, protocol validation, and physical accuracy claims
separate.

1. Contract tests check derived capability tiers, atomic units, electron parity,
   evidence completeness, real-source traceability, and provenance fingerprinting.
2. Replay tests check byte determinism, exact component/operator/overlap round trips,
   coordinate misses, cross-dataset tokens, and independent manifest, array, and
   overlap corruption.
3. Admission tests rerun the v0.22.1 physical SOC symmetry, component-derivative, and
   cross-geometry differential audits through the file-backed provider.
4. Negative controls independently remove reference, basis, method, frame-invariance,
   and tracking evidence; each must fail only the real-admission layer.
5. An unconverged record fails the protocol layer, and a static-only capability fails
   moving-nuclear use.
6. The PySCF boundary fails closed when the runtime or a method-specific validated
   provider is absent.

## Canonical acceptance

`examples/121_recompute_v0230_campaign.py` writes
`results/v0230_molecular_soc_admission_campaign.json` and both reference replay
directories. The canonical result is **93/93 gates passing**:

- 67 inherited v0.22.1 gates;
- 26 new v0.23.0 gates.

The new gates cover capability tiers, both parity sectors, exact replay, deterministic
bytes and fingerprints, integrity corruption, exact-coordinate policy, convergence,
all five evidence families, fixture-versus-real separation, and the PySCF runtime
boundary.

## Reproducibility

Run:

```bash
python -m pytest -q
python examples/121_recompute_v0230_campaign.py
```

The expected dataset fingerprints are:

```text
even: a0c90420ace96c899b5e033a8b03d43cd316eb17511985140541be4e4dec8255
odd:  9e12baf2950ee2f7fe6ec5debe183cc3c5b2ff1e25ef3df20fe636e05a2c9acb
```

## Interpretation

A green v0.23.0 campaign proves that the admission machinery distinguishes valid
trajectory data from missing physical evidence and that captured data replay exactly.
It does not prove that an ab-initio SOC method is accurate. The canonical fixtures are
analytic, and the release environment has no PySCF installation. Consequently:

- molecular SOC protocol validated: **yes**;
- deterministic replay validated: **yes**;
- real molecular SOC backend admitted: **no**;
- ab-initio SOC validated: **no**;
- live PySCF SOC runtime validated: **no**.
