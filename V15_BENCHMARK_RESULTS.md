# v0.15 Benchmark Results

Canonical machine-readable output:

```text
results/v015_cost_aware_cache_campaign.json
```

## 1. Physical result

```text
initial basis size: 10
final basis size: 11
average basis size: 10.925

projection fidelity: 0.8822514544600691
relative initial residual: 0.11774854553993085
initial reduced-density error: 0.033619920355630904

projected-state dynamics density error:
9.527804623556872e-05

original-target density error:
0.03330494031478426

trace distance:
0.02355014914359723

population L2 error:
0.028084897912094255

Gaussian populations:
[0.2458651322299296, 0.7541348677700704]

exact target populations:
[0.22600611046735578, 0.7739938895326441]

Gaussian purity:
0.6613299727370956

exact target purity:
0.6762081969769371

coherence phase error / rad:
0.0028906431794135244

maximum generalized norm drift:
2.115581485107043e-06

maximum condition number:
1470.755892050532
```

## 2. Cost-aware adaptive event

```text
step: 10
time: 0.05

candidate:
parent=9;target=0;pos=target_force:0.06;mom=nac;width=1.35

basis:
10 -> 11

relative defect:
0.03238019095782259
->
0.029987957971150326

predicted capture fraction:
0.1506769024198026

residual-only best:
parent=9;target=0;pos=target_force:0.06;mom=nac;width=1.35

residual-only best capture:
0.1506769024198026

cost-aware utility:
0.24210065440802597

normalized incremental cost:
0.6223729662698816

estimated incremental horizon seconds:
0.03330978168687729

expanded condition number:
1454.4566005006461

candidate count:
560

residual shortlist:
8

new pair factorizations during matrix expansion:
0
```

For this benchmark the residual-only and cost-aware selectors choose the same physical
candidate. That is useful because it isolates the cache/incremental optimization from a
change in the physical basis-growth decision.

## 3. v0.14 physics comparison

Maximum difference across the stored acceptance metrics:

```text
8.412825991399586e-12
```

Per metric:

```json
{
  "coherence_phase_error": 1.3708652268906718e-15,
  "initial_density_error": 0.0,
  "max_condition_number": 8.412825991399586e-12,
  "max_norm_drift": 2.4424906541753444e-15,
  "projected_dynamics_density_error": 4.24237468071853e-15,
  "projection_fidelity": 0.0,
  "purity": 3.3306690738754696e-16,
  "target_density_error": 7.917277944358148e-15,
  "target_population_error": 4.4374226515486725e-15
}
```

The differences are floating-point level.

## 4. Pair-cache complexity

```text
pair requests:
23826

all pair factorizations:
15763

propagation pair factorizations:
15675

candidate-conditioning pair factorizations:
88

direct cache hits:
8063

reverse views:
7205

inherited pairs reused:
440

cache reuse fraction:
0.6408125577100646
```

## 5. v0.14 factorization-equivalent comparison

```text
v0.14 propagation factorization-equivalent baseline:
103103

v0.15 propagation pair factorizations:
15675

avoided:
87428

reduction:
84.797 %
```

This is the central portable performance result.

## 6. Wall-time diagnostic

```text
saved v0.14 adaptive runtime:
11.289004 s

v0.15 adaptive runtime:
4.207306 s

diagnostic speedup:
2.683 x

runtime reduction:
62.73 %
```

Wall time is environment dependent and is not part of release acceptance.

## 7. Timing categories

```text
cached endpoint S/H:
1.133081 s

cached moving-basis T:
1.668603 s

TDSE defect:
0.254956 s

candidate residual ranking:
0.156593 s

cost reranking:
0.000133 s

Cayley solves:
0.006800 s

pruning audits:
0.002573 s

total:
4.207306 s
```

## 8. Acceptance

```json
{
  "checks": {
    "cache_reuse": true,
    "coherence_phase": true,
    "conditioning": true,
    "cost_aware_enrichment": true,
    "cost_utility_gate": true,
    "enrichment_reduces_defect": true,
    "incremental_expansion_reuses_candidate_pairs": true,
    "initial_density_representation": true,
    "norm": true,
    "pair_factorization_reduction": true,
    "projected_dynamics": true,
    "target_density": true,
    "target_population": true,
    "v14_physics_regression": true
  },
  "passed": true,
  "thresholds": {
    "max_coherence_phase_error": 0.0035,
    "max_condition_number": 5000.0,
    "max_incremental_expansion_pair_factorizations": 0,
    "max_initial_density_error": 0.035,
    "max_norm_drift": 0.0001,
    "max_projected_dynamics_density_error": 0.003,
    "max_target_density_error": 0.035,
    "max_target_population_error": 0.03,
    "max_v14_reference_metric_difference": 1e-09,
    "min_cache_hit_fraction": 0.6,
    "min_cost_aware_utility": 0.15,
    "min_enrichment_events": 1,
    "min_factorization_reduction": 0.84
  },
  "v14_reference_difference": {
    "maximum": 8.412825991399586e-12,
    "per_metric": {
      "coherence_phase_error": 1.3708652268906718e-15,
      "initial_density_error": 0.0,
      "max_condition_number": 8.412825991399586e-12,
      "max_norm_drift": 2.4424906541753444e-15,
      "projected_dynamics_density_error": 4.24237468071853e-15,
      "projection_fidelity": 0.0,
      "purity": 3.3306690738754696e-16,
      "target_density_error": 7.917277944358148e-15,
      "target_population_error": 4.4374226515486725e-15
    }
  }
}
```

All configured v0.15 acceptance checks pass.
