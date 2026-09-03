# v0.16 Benchmark Results

Canonical machine-readable output:

```text
results/v016_sparse_locality_campaign.json
```

## Primary physical result

```text
initial basis size:
10

final basis size:
11

average basis size:
10.925

projection fidelity:
0.8822514544600691

initial density error:
0.033619920355630904

projected-state dynamics density error:
0.00013361460054812487

target density error:
0.03333954068459046

target population error:
0.028199413658981914

trace distance:
0.023574615299718705

purity:
0.6612060347623693

coherence phase error:
0.0029095064228781115

maximum norm drift:
2.0053154308197207e-06

maximum condition number:
1431.0606683729504
```

## Sparse adaptive event

```text
step:
10

time:
0.05

candidate:
parent=9;target=0;pos=target_force:0.06;mom=nac;width=1.35

basis:
10 -> 11

relative defect:
0.03255613084426844
->
0.030270053821475414

capture fraction:
0.14896413884583096

local utility:
0.5112460846439298

normalized incremental cost:
0.2913746301833896

predicted local degree:
9

electronic cache hit:
True

active edges after insertion:
52

edge fraction after insertion:
0.9454545454545454
```

## Comparison with dense v0.15

```text
final reduced-density difference:
0.0001631265880682423

v0.16 propagation pair factorizations:
14973

v0.15 propagation pair factorizations:
15675

pair reduction:
4.48 %
```

The compact benchmark remains highly connected, so the pair reduction is intentionally
small.

## Final dense matrix audit

```text
relative S error:
0.005191742661052565

relative H error:
0.003962632349871911

relative nuclear S error:
0.005191742661052565

omitted off-diagonal pairs:
3

maximum omitted overlap:
0.02057350476995086

maximum omitted H block:
0.2318916964441307
```

## Sparse scaling benchmark

```text
N=20
active edges=37
edge fraction=0.194737
pair reduction=72.86%
dense assembly=0.028987 s
sparse assembly=0.004917 s
speedup=5.90 x
```
```text
N=40
active edges=77
edge fraction=0.098718
pair reduction=85.73%
dense assembly=0.112917 s
sparse assembly=0.008293 s
speedup=13.62 x
```
```text
N=80
active edges=157
edge fraction=0.049684
pair reduction=92.69%
dense assembly=0.523068 s
sparse assembly=0.016066 s
speedup=32.56 x
```

At $N=80$:

```text
pair reduction:
92.69 %

H matrix density:
0.061484375

dense assembly:
0.523068 s

sparse assembly:
0.016066 s

diagnostic assembly speedup:
32.56 x
```

Fitted exponents:

```text
active edge exponent:
1.0425836916313385

KD-tree spatial-candidate exponent:
1.0425836916313385

pair-factorization exponent:
1.027926617366759

dense canonical pair exponent:
1.9737662900529327
```

## Electronic cost demonstration

```text
cached geometry:
    cache hit = True
    normalized cost = 1.675

new geometry:
    cache hit = False
    normalized cost = 3.625
```

This demonstrates that a future ab-initio candidate can be penalized for requiring an
uncached electronic-structure point.

## Acceptance

```json
{
  "checks": {
    "coherence_phase": true,
    "conditioning": true,
    "dense_pair_scaling": true,
    "electronic_cache_cost_demo": true,
    "final_sparse_H_audit": true,
    "final_sparse_S_audit": true,
    "graph_is_actually_sparse": true,
    "initial_density_representation": true,
    "local_edge_scaling": true,
    "n80_edge_fraction": true,
    "n80_pair_reduction": true,
    "norm": true,
    "pair_work_reduced_vs_v15": true,
    "projected_dynamics": true,
    "sparse_result_close_to_v15": true,
    "target_density": true,
    "target_population": true
  },
  "final_rho_difference_vs_v15": 0.0001631265880682423,
  "pair_factorization_reduction_vs_v15": 0.04478468899521526,
  "passed": true,
  "thresholds": {
    "electronic_cache_demo_requires_lower_cost": true,
    "max_coherence_phase_error": 0.0035,
    "max_condition_number": 5000.0,
    "max_final_rho_difference_vs_v15": 0.0015,
    "max_final_sparse_H_relative_error": 0.01,
    "max_final_sparse_S_relative_error": 0.01,
    "max_initial_density_error": 0.035,
    "max_local_edge_scaling_exponent": 1.2,
    "max_norm_drift": 0.0001,
    "max_projected_dynamics_density_error": 0.001,
    "max_scaling_edge_fraction_n80": 0.08,
    "max_target_density_error": 0.035,
    "max_target_population_error": 0.03,
    "min_average_graph_sparsity": 0.04,
    "min_dense_pair_scaling_exponent": 1.8,
    "min_pair_factorization_reduction_vs_v15": 0.04,
    "min_scaling_pair_reduction_n80": 0.9
  }
}
```

All configured v0.16 release criteria pass.
