# v0.11 Benchmark Results

This file records the deterministic v0.11 strong-conical-intersection
basis-completeness campaign.

## 1. Exact reference

The primary electronic observable is the reduced electronic density matrix in the
analytic model's global diabatic basis.

```text
exact diabatic populations: [0.22600611046735578, 0.7739938895326441]
exact reduced-state purity: 0.676208196977
exact linear entropy: 0.323791803023
```

## 2. v0.10 baseline

```text
diabatic populations: [0.9502784027695176, 0.04972159723048253]
population L2 error: 1.02427569862
full-density Frobenius error: 1.02934312365
purity: 0.974885247078
purity error: 0.298677050101
maximum norm drift: 8.055582e-08
maximum condition number: 1.992727e+04
```

This is the basis-completeness failure that motivated v0.11.

## 3. v0.11 reference calculation

The 10-TBF v0.11 reference gives

```text
diabatic populations: [0.23511178903234, 0.7648882109676599]
population L2 error: 0.0128773741212
full-density Frobenius error: 0.15991833275
purity: 0.643115591123
purity error: 0.0330926058538
linear entropy: 0.356884408877
maximum norm drift: 5.615808e-03
maximum condition number: 1.550639e+04
```

Relative to the v0.10 baseline:

```text
population-error improvement: 79.54 x
full-density-error improvement: 6.44 x
purity-error improvement: 9.03 x
```

This is a large improvement, but it is not full convergence.

## 4. Basis ladder

| Nmax | final TBFs | population error | density error | purity | max norm drift | max condition |
|---:|---:|---:|---:|---:|---:|---:|
| 2 | 2 | 1.0555481 | 1.0578969 | 0.99990798 | 2.024e-08 | 1.049e+00 |
| 4 | 4 | 1.0565289 | 1.0588022 | 0.99996352 | 3.353e-07 | 8.712e+01 |
| 6 | 6 | 1.0544716 | 1.0570544 | 0.99996 | 6.297e-06 | 1.554e+03 |
| 8 | 8 | 0.96436258 | 0.96567314 | 0.87736867 | 1.871e-05 | 1.417e+04 |
| 10 | 10 | 0.012877374 | 0.15991833 | 0.64311559 | 5.616e-03 | 1.551e+04 |

The convergence is strongly non-monotonic at small basis size. Substantial transfer
appears only once the basis becomes sufficiently rich. Therefore the 10-TBF result
should not be described as a fully basis-converged limit.

## 5. Basis-completeness diagnostics

```text
lineage depth: 3
generation histogram: {'0': 1, '1': 3, '2': 4, '3': 2}
width diversity ratio: 31.855414
canonical participation ratio: 3.3241503
```

The raw basis contains 10 TBFs, but the coefficient-weighted canonical participation
ratio is only about 3.32. Raw TBF count alone
is therefore not used as the definition of basis completeness.

## 6. Ablations

| ablation | population error | density error | purity | max norm drift | max condition |
|---|---:|---:|---:|---:|---:|
| fixed_width_only | 0.58884686 | 0.72614997 | 0.78452499 | 3.267e-01 | 8.987e+06 |
| no_position_optimization | 1.0396142 | 1.0456189 | 0.99996469 | 2.406e-08 | 1.738e+05 |
| single_child_per_event | 0.028785258 | 0.16062125 | 0.71036087 | 1.117e-02 | 1.528e+04 |

The no-position-search ablation returns to an almost pure electronic state and a large
population error. The fixed-width-only ablation is both much less accurate and poorly
conditioned. One child per event already recovers good diagonal populations, but its
full-density error remains of the same order as the two-child reference.

## 7. Acceptance result

```json
{
  "checks": {
    "conditioning": true,
    "full_density": false,
    "norm": true,
    "population": true,
    "purity": true
  },
  "passed": false,
  "thresholds": {
    "max_condition_number": 1000000.0,
    "max_density_frobenius_error": 0.1,
    "max_norm_drift": 0.01,
    "max_population_l2_error": 0.05,
    "max_purity_error": 0.05
  }
}
```

The v0.11 reference passes the diagonal-population, reduced-purity, generalized-norm,
and overlap-conditioning checks. It fails the configured full reduced-density
criterion.

**Conclusion:** v0.11 achieves partial basis convergence, not full quantum
convergence. The remaining discrepancy is concentrated more strongly in electronic
coherence/phase structure than in the diagonal populations.

## 8. Reproducibility

The complete machine-readable campaign is:

```text
results/v011_basis_completeness_campaign.json
```

The expensive ablations were executed in isolated processes and merged into this
file. This avoids retaining multiple large dynamic electronic gauge graphs in memory
simultaneously and does not alter the individual calculations.
