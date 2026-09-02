# v0.13 Build Validation Report

Validated on 2026-08-12.

## Source validation

All 193 Python files in the repository parsed successfully with Python's AST parser.

## Automated regression suite

```text
158 passed in 8.19 s
```

The cumulative suite includes all retained v0.1-v0.12 tests plus new v0.13 checks for:

- monotonic Hilbert residual reduction;
- exact one-candidate predicted residual gain versus actual reprojection gain;
- deterministic density-screened residual selection;
- direct and vectorized/prepared residual builders selecting the same pure-greedy
  candidates on a controlled dictionary;
- finite instantaneous TDSE defect reconstruction;
- near-orthogonality of the TDSE defect to the current represented basis;
- positive residual capture by an admissible defect candidate;
- zero-coefficient insertion preserving the instantaneous wavefunction;
- actual defect reduction after Galerkin-space enlargement;
- predicted versus actual squared-defect reduction;
- acceptance logic that independently checks residual monotonicity and defect-gain
  prediction.

## Release campaign

The complete v0.13 release campaign was executed successfully and saved as:

```text
results/v013_residual_driven_campaign.json
```

The candidate dictionary contains:

```text
726 normalized Gaussian candidates
position radius = 1.0
position spacing = 0.2
width scales = [1.0, 1.5, 2.0, 3.0, 4.0, 6.0]
```

## v0.13 reference result

```text
basis size: 11
projection fidelity: 0.8902521956060818
relative residual: 0.10974780439391818
initial density error: 0.03209140317550961

projected-state dynamics density error:
0.00011354880287339317

original-target density error:
0.03178630139393256

target population error:
0.025521902605714804

trace distance:
0.022476309264489125

purity error:
0.0126522377254491

coherence phase error / rad:
0.0023799927838891125

maximum generalized norm drift:
1.0593904686828637e-06

maximum overlap condition number:
3465.8914579773386
```

## v0.12 comparison

```text
v0.12 projection fidelity:        0.832276023595292
v0.13 projection fidelity:        0.8902521956060818

v0.12 relative residual:          0.16772397640470793
v0.13 relative residual:          0.10974780439391818

v0.12 target density error:       0.03500028070905269
v0.13 target density error:       0.03178630139393256

v0.12 projected dynamics error:   0.00029022869338069174
v0.13 projected dynamics error:   0.00011354880287339317
```

The v0.13 basis is not manually placed. It is selected from a documented deterministic
dictionary by residual reduction followed by an initial-density screen restricted to
the top residual candidates.

## TDSE-defect enrichment

```text
selected candidate:
dq=(0.0, -0.4);dp=(0.0, 0.0);width_scale=4

defect norm before:
0.31502411763651

defect norm after:
0.28652498129723825

predicted squared reduction:
0.01714362854505765

actual squared reduction:
0.017143629785278974

capture fraction:
0.17274884030759896

expanded condition number:
7460.863531968121
```

The relative error between predicted and actual squared-defect reduction is:

```text
7.234298863186708e-08
```

The enrichment candidate is inserted with a zero two-component electronic coefficient,
so the physical wavefunction is unchanged at the insertion instant.

## Acceptance result

```json
{
  "checks": {
    "coherence_phase": true,
    "conditioning": true,
    "defect_gain_prediction": true,
    "defect_reduction": true,
    "initial_density_representation": true,
    "monotone_residual_refinement": true,
    "norm": true,
    "projected_dynamics": true,
    "target_full_density": true,
    "target_population": true
  },
  "defect_prediction_relative_error": 7.234298863186708e-08,
  "passed": true,
  "thresholds": {
    "max_coherence_phase_error": 0.003,
    "max_condition_number": 5000.0,
    "max_defect_prediction_relative_error": 0.005,
    "max_initial_density_error": 0.033,
    "max_norm_drift": 0.0001,
    "max_projected_dynamics_density_error": 0.0002,
    "max_target_density_error": 0.033,
    "max_target_population_error": 0.03,
    "min_defect_squared_reduction": 1e-08
  }
}
```

All configured v0.13 release criteria pass.

## Representative examples validated

The following release-facing examples were executed successfully:

```text
examples/48_v013_residual_selection.py
examples/49_v013_projection_ladder.py
examples/50_v013_tdse_defect.py
examples/51_v012_v013_comparison.py
```

`examples/52_recompute_v013_campaign.py` contains the full campaign recomputation
workflow. The canonical release campaign itself was executed during the build and
stored in `results/v013_residual_driven_campaign.json`.

## PySCF status

The explicit PySCF SA-CASSCF backend, many-electron cross-geometry overlaps, state
tracking, and gauge graph from earlier releases remain in the repository and in the
cumulative regression suite.

The v0.13 residual benchmark uses the analytic 2D LVC model. A full-dimensional
molecular TDSE residual is not claimed.

See:

```text
V13_PYSCF_RESIDUAL_BRIDGE.md
```

for the proposed sampled local-diabatic molecular extension.
