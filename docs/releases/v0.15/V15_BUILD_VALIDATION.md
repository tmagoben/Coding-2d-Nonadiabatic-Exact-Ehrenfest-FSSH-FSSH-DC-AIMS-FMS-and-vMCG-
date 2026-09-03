# v0.15 Build Validation Report

Validated on 2026-08-12.

## Source validation

All 226 Python files parsed successfully with Python's AST parser before packaging.

## Automated regression suite

```text
178 passed in 9.10 s
```

The suite retains all earlier release tests and adds v0.15 checks for:

- cached pair algebra against the original analytic S/H/T implementations;
- exact reverse-orientation pair identities;
- no extra pair factorization when T reuses an already primed S/H cache;
- incremental accepted-child matrix expansion against a full rebuild;
- incremental pair-cache subset/pruning against a full rebuild;
- explicit v0.14/v0.15 factorization-equivalent counts;
- cost-model pair and dense-solve growth;
- condition-aware cost reranking;
- cost-utility rejection of low-benefit candidates;
- fixed-basis v0.15 propagation against v0.14;
- real short-run cost-aware TDSE-defect enrichment;
- v0.15 release acceptance logic.

## Canonical release campaign

```text
results/v015_cost_aware_cache_campaign.json
```

All configured acceptance checks pass.

## Physical invariance

Maximum difference between the stored v0.15 and v0.14 reference acceptance metrics:

```text
8.412825991399586e-12
```

The release threshold is `1e-9`.

The final target-density errors are:

```text
v0.14: 0.03330494031479218
v0.15: 0.03330494031478426
```

The cache/incremental optimization therefore preserves the reference dynamics to
floating-point precision.

## Adaptive event

```text
step: 10
time: 0.05

candidate:
parent=9;target=0;pos=target_force:0.06;mom=nac;width=1.35

basis:
10 -> 11

relative TDSE defect:
0.03238019095782259
->
0.029987957971150326

capture fraction:
0.1506769024198026

cost-aware utility:
0.24210065440802597

normalized incremental cost:
0.6223729662698816

estimated incremental horizon seconds:
0.03330978168687729

new pair factorizations during matrix expansion:
0
```

The candidate enters with zero electronic amplitude.

## Pair-cache audit

```text
pair requests:
23826

all pair factorizations:
15763

propagation pair factorizations:
15675

candidate-conditioning pair factorizations:
88

direct hits:
8063

reverse views:
7205

inherited pair reuse:
440

cache reuse fraction:
0.6408125577100646
```

## v0.14 factorization-equivalent comparison

```text
v0.14 propagation baseline:
103103

v0.15 propagation factorizations:
15675

avoided:
87428

reduction:
84.797 %
```

The release requirement is at least 84%.

## Timing diagnostic

Saved benchmark timings:

```text
v0.14 adaptive runtime:
11.289004 s

v0.15 adaptive runtime:
4.207306 s

diagnostic speedup:
2.683 x

diagnostic runtime reduction:
62.73 %
```

Wall time is not an acceptance criterion.

Measured v0.15 timing categories:

```text
cached endpoint S/H:
1.133081 s

cached moving-basis T:
1.668603 s

TDSE defect:
0.254956 s

candidate residual ranking:
0.156593 s

cost-aware reranking:
0.000133 s

Cayley solves:
0.006800 s

total adaptive run:
4.207306 s
```

## Reference accuracy

```text
initial reduced-density error:
0.033619920355630904

projected-state dynamics density error:
9.527804623556872e-05

original-target density error:
0.03330494031478426

population L2 error:
0.028084897912094255

coherence phase error / rad:
0.0028906431794135244

maximum generalized norm drift:
2.115581485107043e-06

maximum condition number:
1470.755892050532
```

## Representative examples executed

```text
examples/59_v015_pair_cache.py
examples/60_v015_cost_aware_event.py
examples/61_v015_complexity.py
examples/62_v014_v015_comparison.py
examples/63_v015_incremental_matrix.py
```

`examples/64_recompute_v015_campaign.py` is the full user-facing recomputation wrapper.

## PySCF status

The inherited PySCF SA-CASSCF, many-electron overlap, state-tracking, and gauge-graph
layers remain in the repository and cumulative regression suite.

v0.15 optimizes the analytic Gaussian dynamics layer. It does not yet include an
electronic-structure-aware molecular cost model.

See `V15_PYSCF_COST_BRIDGE.md`.
