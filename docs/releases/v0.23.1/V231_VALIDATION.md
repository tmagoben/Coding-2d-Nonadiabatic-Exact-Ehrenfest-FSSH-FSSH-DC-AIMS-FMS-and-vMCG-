# v0.23.1 validation

## Canonical campaign

`examples/122_recompute_v0231_campaign.py` builds both parity-sector dossiers and
writes `results/v0231_raw_evidence_admission_campaign.json`. The campaign passes
**123/123 gates**: 93 inherited v0.23.0 gates and 30 new gates.

The new gates cover:

- deterministic dossier bytes, fingerprints, and artifact inventories;
- even singlet–triplet and odd two-doublet protocol paths;
- exact replay fingerprint and receipt-coordinate binding;
- complete per-stage convergence;
- derived reference, basis, method, frame, and connected-manifold tracking evidence;
- exact agreement between raw-derived evidence and provider provenance;
- raw artifact and dossier corruption;
- missing, unconverged, duplicate-output, wrong-coordinate, and wrong-environment
  receipts;
- tampered evidence summaries, disconnected tracking graphs, and path traversal;
- synthetic external relabelling that v0.23.0 summary evidence alone would accept;
- mandatory backend parser attestation and executable validation;
- mandatory fresh runtime execution for live admission;
- PySCF unavailable/incomplete fail-closed behavior.

## Canonical fingerprints

```text
even replay   2cabad829ff634ca5612c1e9e8acfbec0afe0a3f6c473ff951fe1fcf73daa651
even dossier  f7d1c3479911ac2b864a12701397eeb99100f0e2038eb4b5a2d0b69b4ffbc82d
odd replay    7cef09125fd86aa0970ee965dcf0ed19448d1d4b91c3fd86341f390d09086797
odd dossier   865d7e88d62b937bef9aee8faf0ffee64031717524d431a3ccef9c561b612487
```

## Interpretation

The campaign demonstrates that raw artifacts and receipts deterministically reproduce
their evidence summaries and that a relabelled fixture cannot cross the new executable
validator boundary. It does not demonstrate molecular accuracy.

- raw-evidence admission protocol validated: **yes**;
- calculation-receipt integrity validated: **yes**;
- external molecular SOC snapshot admitted: **no**;
- live molecular SOC backend admitted: **no**;
- ab-initio SOC validated: **no**;
- live PySCF SOC runtime validated: **no**.
