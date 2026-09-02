# v0.21.4 Validation Contract

Machine-readable campaign:

```text
results/v0214_differential_restart_campaign.json
```

All **21/21** v0.21.4 gates pass, including inherited v0.21.3 acceptance. The cumulative
source suite contains **280 passing tests**.

## Differential provider certification

| Fixture | K scaled error | D scaled error | Max overlap-isometry residual |
|---|---:|---:|---:|
| fixed frame | 8.232012296296499e-15 | 0.0 | 0.0 |
| coordinate-dependent complex frame | 6.566733600790935e-14 | 1.6631553311401476e-11 | 2.4405271014360246e-16 |

The fixed-frame and moving-frame providers pass the structural, K, D, overlap, and
provenance checks. Two negative controls remain pointwise structurally valid: adding a
Hermitian $10^{-3}I$ defect to K fails only the physical-derivative gate, while erasing
the nonzero connection fails only the D gate.

## Zero-SOC rehearsal

H, K, D, mass, and cross-geometry overlap errors are exactly zero across the campaign
geometries. The rehearsal provider also passes the differential contract. Spin-free and
explicit-zero-SOC trajectories give

| Quantity | Error |
|---|---:|
| final positions | 0.0 |
| final momenta | 0.0 |
| phase-aligned metric coefficient vector | 7.255325009694841e-17 |

## Checkpoint/restart

| Path | Position error | Momentum error | Coefficient error |
|---|---:|---:|---:|
| dense fixed frame | 0.0 | 0.0 | 8.238824461827495e-16 |
| sparse fixed frame | 0.0 | 0.0 | 8.238824461827495e-16 |
| dense moving complex frame | 0.0 | 0.0 | 1.9932132794842223e-15 |

The saved digest recomputes exactly after round trip. A one-element coordinate mutation
without a matching digest is rejected. Sparse edge `(3, 8)` survives the segment
boundary, and uninterrupted/resumed final edge sets agree.

An exactly zero local coefficient block retains its valid guide density across restart.
An adaptive run inserts UID 13 at global step 2, checkpoints it and its inherited guide
density, then prunes it at global step 3 after resume; the final live UIDs are `[3, 8]`.

## Acceptance gates

The 21 gates cover:

1. fixed-frame differential consistency;
2. coordinate-dependent complex-frame differential consistency;
3. wrong-K negative-control detection;
4. wrong-D negative-control detection;
5. exact zero-SOC H equivalence;
6. exact zero-SOC K equivalence;
7. exact zero-SOC D/mass/overlap equivalence;
8. zero-SOC differential consistency;
9. zero-SOC position equivalence;
10. zero-SOC momentum equivalence;
11. zero-SOC coefficient equivalence;
12. dense restart positions;
13. dense restart momenta;
14. dense restart coefficients;
15. checkpoint digest round trip;
16. checkpoint corruption rejection;
17. sparse graph restart;
18. moving complex-frame restart;
19. zero-block guide-density restart;
20. global adaptive lifecycle restart;
21. inherited v0.21.3 acceptance.

## Scope boundary

The campaign contains no physical H_SOC or K_SOC and makes no real PySCF runtime claim.
Passing establishes analytic interface and deterministic-restart readiness. It does not
establish SOC dynamics, ab-initio SOC accuracy, exact-grid SOC agreement, production
AIMS, or production asynchronous scheduling.
