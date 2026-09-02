# v0.22.1 validation record

## Acceptance structure

The canonical campaign contains **67 gates**:

- 53 inherited v0.22.0 gates, unchanged;
- 14 new derivative, contract, symmetry, exact-grid, convergence, and serialization
  gates.

All **67/67** gates pass. The machine-readable record is
`results/v0221_corrective_hardening_campaign.json`.

## Full-matrix derivative audit

The arbitrary three-state, two-coordinate provider passes without a `.config` object:

| Quantity | Maximum residual |
|---|---:|
| Spin-free component derivative | 4.0135369935227635e-14 |
| SOC component derivative | 2.51964190546808e-15 |

Both are below the `2e-9` acceptance threshold. Six rows are audited: three finite-
difference steps for each of two nuclear coordinates.

The cancelling-error fixture preserves H/K composition, the total transported
differential, and the sampled scalar SOC-force check, but is correctly rejected by the
separate component checks. Its spin-free and SOC residuals are respectively
`6.686815952118713e-4` and `6.686815952118393e-4`.

## Symmetry negative controls

- A mixed singlet/doublet model space is rejected by electron-parity consistency.
- A nonunitary time-reversal matrix with exact `JJ* = I` is rejected; its unitarity
  residual is `0.8401680504168058` while its square residual is zero.
- Swapped numerical projectors in provenance are rejected with symmetry-provenance
  residual `1.1547005383792517`.
- The unchanged singlet–triplet reference passes every symmetry-admission condition.

## Exact-grid corrections

A five-step run with `store_every=2` records times

```text
[0.00, 0.02, 0.04, 0.05]
```

so the true final time is retained. A wrapper with no `.config` passes using the mass
emitted by the electronic contract. Fixed frame and constant scalar mass are certified;
a coordinate-dependent electronic frame is rejected. The optimized precomputed
trajectory is also compared directly with repeated calls to the public step operator.

All inherited v0.22.0 exact-grid norm, energy, timestep, spatial-resolution, box, and
Gaussian/grid population gates continue to pass.

## Gaussian-basis SOC convergence

At `t=1.0`, the physical triplet population is:

| Gaussian basis size | Initial projection fidelity | Final triplet population |
|---:|---:|---:|
| 1 | 0.9999999999999998 | 1.1589943670061445e-05 |
| 3 | 0.9999999999999998 | 1.1774963819908163e-05 |
| 5 | 0.9999999999999998 | 1.1775857039484610e-05 |

The 1-to-3 difference is `1.8502014984671707e-7`; the 3-to-5 difference is
`8.932195764469294e-10`; the narrowing ratio is `0.004827688104170986`.

## Sparse-threshold SOC convergence

| Entry threshold | Active edges | Metric coefficient error vs dense |
|---:|---:|---:|
| 1.20 | 0 | 1.4803870179355646e-03 |
| 0.50 | 9 | 7.248410605232162e-04 |
| 0.15 | 14 | 2.0668976179392948e-04 |
| 0.05 | 15 | 8.673617405127575e-19 |

The edge count grows monotonically and the trajectory error falls monotonically to
the dense result.

## Claim boundary

These results validate the implemented physical analytic SOC models and framework
contracts. They do not validate molecular or ab-initio SOC matrix elements, a PySCF
SOC runtime, external magnetic fields, production AIMS equations, or molecular
predictive accuracy.

