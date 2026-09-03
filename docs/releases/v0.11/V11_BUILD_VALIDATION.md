# v0.11 Build Validation Report

Validated on 2026-08-12.

## Source validation

All 155 repository Python files parse successfully with Python's AST parser.

## Automated regression suite

```text
[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m [ 57%]
[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m                    [100%][0m
[32m[32m[1m125 passed[0m[32m in 5.42s[0m[0m
```

The cumulative suite includes every retained v0.1-v0.10 regression plus new v0.11
checks for:

- exact unequal-width multidimensional Gaussian overlap;
- unequal-width gradient and kinetic matrix elements against 2D quadrature;
- unequal-width moving-basis derivatives, including an explicit width derivative,
  against finite differences;
- unequal-width graph-SPA Hermiticity;
- energy-conserving shifted child positions and momenta;
- deterministic spawn candidate ranking;
- sibling and existing-basis redundancy blocking;
- multi-generation v0.11 propagation;
- reduced electronic density with unequal widths;
- canonical nonorthogonal-basis participation diagnostics;
- release acceptance that distinguishes populations from the full reduced density;
- campaign JSON serialization of complex density matrices.

## Representative executable outputs

### Unequal-width Gaussian algebra

```text
Unequal-width Gaussian algebra
------------------------------
<g_i|g_j> = (0.44275496381026447-0.036969415783674964j)
real overlap saddle = [-0.15264765 -0.18523529]
complex cross centroid = [-0.15264765-0.45425834j -0.18523529+0.4258343j ]
<g_i|T|g_j> = (0.05816007681977689-0.007663499691758646j)
```

### Optimized-spawning-inspired candidate search

```text
Top optimized-spawn-inspired local candidates
---------------------------------------------
 1 score=7.276851e-03 coupling=7.276851e-03 shift=-0.050(nac) width_scale=0.650 momentum=nac |S_nuc|=0.792367 dE=-2.082e-17
 2 score=7.154505e-03 coupling=7.154505e-03 shift=+0.050(momentum) width_scale=0.650 momentum=nac |S_nuc|=0.802074 dE=-6.939e-18
 3 score=7.097003e-03 coupling=7.097003e-03 shift=+0.050(target_force) width_scale=0.650 momentum=nac |S_nuc|=0.806222 dE=-6.939e-18
 4 score=7.009861e-03 coupling=7.009861e-03 shift=+0.000(none) width_scale=0.650 momentum=nac |S_nuc|=0.813009 dE=-6.939e-18
 5 score=6.909765e-03 coupling=6.909765e-03 shift=-0.050(target_force) width_scale=0.650 momentum=nac |S_nuc|=0.818861 dE=-6.939e-18
 6 score=6.839482e-03 coupling=6.839482e-03 shift=-0.050(momentum) width_scale=0.650 momentum=nac |S_nuc|=0.823295 dE=-2.082e-17
 7 score=6.727670e-03 coupling=6.727670e-03 shift=+0.050(nac) width_scale=0.650 momentum=nac |S_nuc|=0.829771 dE=-6.939e-18
 8 score=6.373012e-03 coupling=6.373012e-03 shift=-0.050(nac) width_scale=1.000 momentum=nac |S_nuc|=0.840907 dE=-2.082e-17
```

### Strong-CI reference run

This executable was also run successfully during the final validation pass:

```text
final basis size: 10
lineage depth: 3
generation histogram: {0: 1, 1: 3, 2: 4, 3: 2}
width diversity ratio: 31.855414208434606
canonical participation ratio: 3.3241503446209557
diabatic reduced populations: [0.23511179 0.76488821]
reduced-state purity: 0.643115591123172
final generalized norm: 0.9943900700000047
maximum recorded condition number: 15506.39347490685
```

### Saved basis ladder

```text
v0.11 saved basis ladder
------------------------
exact populations: [0.22600611046735578, 0.7739938895326441]
Nmax= 2 Nfinal= 2 population_error=1.055548e+00 density_error=1.057897e+00 purity=0.99990798 norm_drift=2.024e-08 cond=1.049e+00
Nmax= 4 Nfinal= 4 population_error=1.056529e+00 density_error=1.058802e+00 purity=0.99996352 norm_drift=3.353e-07 cond=8.712e+01
Nmax= 6 Nfinal= 6 population_error=1.054472e+00 density_error=1.057054e+00 purity=0.99996000 norm_drift=6.297e-06 cond=1.554e+03
Nmax= 8 Nfinal= 8 population_error=9.643626e-01 density_error=9.656731e-01 purity=0.87736867 norm_drift=1.871e-05 cond=1.417e+04
Nmax=10 Nfinal=10 population_error=1.287737e-02 density_error=1.599183e-01 purity=0.64311559 norm_drift=5.616e-03 cond=1.551e+04

This example reads the release campaign instead of recomputing several large dynamic gauge graphs. Use example 41 to regenerate the base campaign.
```

### Saved ablation study

```text
v0.11 saved branching ablation study
------------------------------------
v0.10 baseline population error: 1.0242756986247679
v0.11 reference population error: 0.012877374121210683
fixed_width_only             population_error=5.888469e-01 density_error=7.261500e-01 purity=0.78452499 norm_drift=3.267e-01 cond=8.987e+06
no_position_optimization     population_error=1.039614e+00 density_error=1.045619e+00 purity=0.99996469 norm_drift=2.406e-08 cond=1.738e+05
single_child_per_event       population_error=2.878526e-02 density_error=1.606213e-01 purity=0.71036087 norm_drift=1.117e-02 cond=1.528e+04

Acceptance: {'checks': {'conditioning': True, 'full_density': False, 'norm': True, 'population': True, 'purity': True}, 'passed': False, 'thresholds': {'max_condition_number': 1000000.0, 'max_density_frobenius_error': 0.1, 'max_norm_drift': 0.01, 'max_population_l2_error': 0.05, 'max_purity_error': 0.05}}
```

## Strong-CI acceptance result

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

The release reference passes the configured diagonal-population, purity, norm, and
conditioning checks but fails the full reduced-density threshold.

That is recorded as **partial convergence**, not full validation.

## Release benchmark files

The primary machine-readable campaign is:

```text
results/v011_basis_completeness_campaign.json
```

The isolated ablation JSON files are retained as provenance for the merged campaign.

## PySCF

The inherited explicit PySCF backend, many-electron SA-CASSCF state tracking, and gauge
graph interfaces remain in the repository. v0.11's release benchmark itself uses the
analytic 2D CI model so the basis-completeness changes can be diagnosed independently
of electronic-structure uncertainty.
