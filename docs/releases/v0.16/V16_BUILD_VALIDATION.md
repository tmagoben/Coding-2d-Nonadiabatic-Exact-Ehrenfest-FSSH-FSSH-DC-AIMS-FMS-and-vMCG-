# v0.16 Build Validation Report

Validated on 2026-08-12.

## Source validation

All 244 Python files in the repository parsed successfully with Python's AST parser.

## Automated regression suite

```text
189 passed in 9.68 s
```

The cumulative suite retains all earlier release tests and adds v0.16 checks for:

- conservative overlap-bound correctness on randomized SPD Gaussian pairs;
- exact far-pair rejection without a Gaussian pair solve;
- persistent enter/exit graph hysteresis;
- full-graph sparse S/H equality with dense cached S/H;
- full-graph sparse T and sparse Cayley equality with the dense moving-basis reference;
- electronic geometry-cache hit/miss cost accounting;
- local-degree cost estimation;
- provider cost changing candidate ordering when appropriate;
- real short-run sparse TDSE-defect enrichment;
- sparse generalized norm conservation;
- release acceptance requiring actual graph sparsity and dense-reference agreement.

## Canonical release campaign

The complete v0.16 release campaign was executed and saved as:

```text
results/v016_sparse_locality_campaign.json
```

All configured release acceptance checks pass.

## Primary sparse-adaptive result

```text
initial basis size:
10

final basis size:
11

average basis size:
10.925

projected-state dynamics density error:
0.00013361460054812487

target density error:
0.03333954068459046

population L2 error:
0.028199413658981914

coherence phase error / rad:
0.0029095064228781115

maximum generalized norm drift:
2.0053154308197207e-06

maximum condition number:
1431.0606683729504
```

## Adaptive event

```text
step: 10
time: 0.05

candidate:
parent=9;target=0;pos=target_force:0.06;mom=nac;width=1.35

relative defect:
0.03255613084426844
->
0.030270053821475414

predicted capture fraction:
0.14896413884583096

local utility:
0.5112460846439298

predicted local degree:
9

active edges after insertion:
52

zero coefficient insertion:
True
```

## Compact-benchmark locality

```text
average off-diagonal graph sparsity:
0.053816793893129766

propagation pair factorizations:
14973

pair-factorization reduction versus v0.15:
0.04478468899521526

final density-matrix difference versus v0.15:
0.0001631265880682423
```

The compact 10–11 TBF CI basis is highly overlapping. The release therefore does not
claim a large sparse speedup on this particular trajectory.

## Dense endpoint sparse-matrix audit

```text
relative S Frobenius error:
0.005191742661052565

relative H Frobenius error:
0.003962632349871911

relative nuclear S Frobenius error:
0.005191742661052565

omitted off-diagonal Gaussian pairs:
3

maximum omitted exact overlap:
0.02057350476995086

maximum omitted H-block Frobenius norm:
0.2318916964441307
```

This audit is intentionally retained because overlap locality does not by itself
constitute a universal Hamiltonian truncation bound.

## Bounded-locality scaling audit

At `N=80`:

```text
active edges:
157

all possible off-diagonal pairs:
3160

KD-tree spatial candidates:
157

globally screened pairs:
3003

pair-factorization reduction:
0.9268518518518518

H matrix density:
0.061484375

dense matrix assembly:
0.523068 s

sparse matrix assembly:
0.016066 s

diagnostic assembly speedup:
32.558 x
```

Fitted exponents over `N=20,40,80`:

```text
active edges:
1.0425836916313385

KD-tree spatial candidates:
1.0425836916313385

pair factorizations:
1.027926617366759

dense canonical pair count:
1.9737662900529327
```

These exponents apply only to the bounded-locality synthetic chain.

## Electronic-structure cost interface

Deterministic geometry-cache demonstration:

```text
cached geometry:
{'cache_hit': True, 'cost_units': 0.05, 'normalized_incremental_cost': 1.675, 'q': [0.04, 0.0]}

new geometry:
{'cache_hit': False, 'cost_units': 2.0, 'normalized_incremental_cost': 3.625, 'q': [1.0, 0.0]}
```

The analytic release itself uses zero provider cost.

## Complexity/timing ledger

```text
endpoint graph updates:
122

midpoint graph updates:
120

exact pair checks:
12395

globally screened pairs:
0

all pair factorizations:
15061

propagation pair factorizations:
14973

candidate pair factorizations:
88

peak active edges:
52

peak S nnz:
230

peak H nnz:
460

graph time:
1.598122 s

sparse S/H assembly:
0.447088 s

sparse T assembly:
0.881072 s

sparse Cayley solves:
0.045959 s

TDSE-defect work:
0.259131 s

total adaptive run:
4.391310 s
```

Wall times are diagnostic and environment dependent.

## Representative examples executed

```text
examples/65_v016_locality_graph.py
examples/66_v016_sparse_audit.py
examples/67_v016_sparse_scaling.py
examples/68_v016_electronic_cost.py
examples/69_v015_v016_comparison.py
```

`examples/70_recompute_v016_campaign.py` contains the complete user-facing campaign
recomputation workflow. The canonical campaign itself was executed directly during the
build.

## Dependencies

v0.16 makes SciPy a core dependency because the release uses:

```text
scipy.spatial.cKDTree
scipy.sparse
scipy.sparse.linalg.spsolve
```

The package metadata now requires:

```text
numpy>=1.24
scipy>=1.10
```

## PySCF status

The inherited PySCF SA-CASSCF, many-electron overlap, state-tracking, and gauge-graph
infrastructure remains in the repository and cumulative regression suite.

v0.16 adds a provider-cost interface but does not claim a calibrated PySCF wall-time
model or a production sparse molecular residual controller.

See `V16_PYSCF_SPARSE_COST_BRIDGE.md`.
