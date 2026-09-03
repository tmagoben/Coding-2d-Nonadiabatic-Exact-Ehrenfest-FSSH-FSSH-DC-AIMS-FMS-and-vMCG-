# v0.23.1 release notes

v0.23.1 hardens molecular spin-orbit-coupling (SOC) admission by replacing
self-contained summary evidence with a traceable raw-evidence dossier. It does not
admit a real backend. Every canonical source remains an explicitly labelled analytic
fixture, and the unavailable PySCF runtime remains a hard stop.

## Added

- Per-calculation receipts binding record role, exact coordinate, backend/method/basis,
  input and output artifacts, and five independent convergence stages.
- SHA-256 and byte-size verification of the calculation template, environment lock,
  every rendered input, every raw output, the independent reference, and any runtime
  probe.
- Reference error, adjacent basis/method changes, translation/rotation residuals, and
  physical-manifold tracking scores derived from stored observations.
- Connected tracking graphs using minimum singular values within complete physical
  manifolds and spectral-norm leakage to competing manifolds.
- Exact binding between the raw dossier, the v0.23.0 replay fingerprint, every replay
  coordinate, and the summarized v0.23.0 evidence embedded in provider provenance.
- Separate external-snapshot and live-backend admission outcomes.
- Executable backend-validator requirement. A dossier attestation alone cannot admit a
  real source; a method-specific parser must run and return exactly `True`.
- A pinned optional PySCF 2.13.1 boundary with the documented
  `state=(ket,bra) -> <bra|d ket/dR>` NAC convention.
- Deterministic even singlet–triplet and odd two-Kramers-doublet raw-evidence bundles.

## Acceptance

The canonical campaign passes **123/123 gates**:

- 93 inherited v0.23.0 gates;
- 30 new receipt, evidence, integrity, tracking, identity, and runtime gates.

Reference fingerprints:

| Sector | Replay fingerprint | Dossier fingerprint |
|---|---|---|
| Even singlet–triplet | `2cabad829ff634ca5612c1e9e8acfbec0afe0a3f6c473ff951fe1fcf73daa651` | `f7d1c3479911ac2b864a12701397eeb99100f0e2038eb4b5a2d0b69b4ffbc82d` |
| Odd two-doublet | `7cef09125fd86aa0970ee965dcf0ed19448d1d4b91c3fd86341f390d09086797` | `865d7e88d62b937bef9aee8faf0ffee64031717524d431a3ccef9c561b612487` |

Recompute with:

```bash
python examples/122_recompute_v0231_campaign.py
```

## Claim boundary

v0.23.1 validates the raw-evidence admission machinery and proves that a synthetic
dataset relabelled as external can no longer pass solely on summary evidence. It does
not validate molecular SOC accuracy, admit an external ab-initio snapshot, admit a live
backend, or validate a live PySCF runtime. Those require real raw outputs and an
executable method-specific validator in the pinned environment.
