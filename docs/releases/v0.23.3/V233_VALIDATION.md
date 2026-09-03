# v0.23.3 validation

The canonical campaign has 208 native Boolean gates:

| Layer | Gates | Purpose |
|---|---:|---|
| Inherited v0.23.2 | 168 | Full prior release, including real spin-free PySCF runtime evidence |
| Finite-manifold transport | 11 | Contraction, polar transport, covariance, retention, rank, reciprocity |
| Replay and NAC compatibility | 15 | Format 2, deterministic migration, quarantine, identity binding |
| Complete manifolds and SOC convention | 10 | Singlet/triplet, doublets, gauges, leakage, Kramers, order/prefactor |
| Runtime profiles | 4 | Canonical lock, compatibility, portability distinction, fingerprints |
| **Total** | **208** | All must pass |

The negative controls intentionally exercise spectral expansion, rank loss,
low retention, broken reciprocity, raw-overlap substitution, missing/mismatched
NAC identity, unknown/wrong-sign legacy data, incomplete projectors, competing
manifold leakage, broken Kramers covariance, SOC state-order/prefactor mismatch,
and portable runtimes that must not claim canonical byte identity.

The campaign is recomputed from source by
`examples/124_recompute_v0233_campaign.py`; its canonical JSON record is stored in
`results/v0233_transport_compatibility_campaign.json`.

Validated claims stop at framework transport/convention machinery, physical
analytic SOC, and the inherited real spin-free PySCF runtime. No external or live
molecular-SOC source is admitted.
