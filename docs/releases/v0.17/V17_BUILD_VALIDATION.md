# v0.17 Build Validation Report

Validated on 2026-08-12.

## Source validation

All 260 Python files in the repository parsed successfully with Python's AST parser.

## Automated regression suite

```text
200 passed in 19.30 s
```

The cumulative suite retains all v0.1-v0.16 tests and adds v0.17 checks for:

- safe global overlap-radius screening;
- pair-specific conservative overlap bounds;
- exact local S/H/T edge-importance construction;
- score hysteresis;
- one-sided score and search-floor relaxation;
- global local-score L2 budget promotion;
- real online dense-audit/rebuild behavior;
- sparse TDSE-defect enrichment under the v0.17 graph;
- threshold-convergence sweeps;
- local-budget-convergence sweeps;
- release acceptance logic.

## Canonical release campaign

```text
results/v017_sparse_error_control_campaign.json
```

All configured release acceptance checks pass.

## Online controller validation

The release intentionally starts with an over-aggressive sparse graph.

Initial audit:

```text
enter score:
0.060

S relative error:
0.01647878134803449

H relative error:
0.016940324970549453

status:
FAIL
```

The controller automatically relaxes to:

```text
enter score:
0.030

exit score:
0.015

search overlap floor:
5e-6
```

Immediate re-audit:

```text
S relative error:
0.0026588687049987378

H relative error:
0.0020460833614795023

status:
PASS
```

All later scheduled audits pass and there are zero unresolved audits.

## Final matrix audit

```text
relative S error:
0.005191742661052565

relative H error:
0.003962632349871911

relative Snuc error:
0.005191742661052565

omitted off-diagonal Gaussian pairs:
3
```

The release tolerance for each matrix error is `0.006`.

## Local importance budget

```text
configured B_local:
<= 0.08

maximum recorded trajectory B_local:
0.024083537794312975
```

The local score budget is a proxy and is not treated as a substitute for the dense
matrix audit.

## Physical result

```text
initial basis size:
10

final basis size:
11

average basis size:
10.925

projected-state dynamics error:
0.00013361460054442858

target density error:
0.03333954068459557

target population error:
0.02819941365898425

coherence phase error:
0.0029095064228609707

maximum generalized norm drift:
2.0053154399235495e-06

maximum condition number:
1431.0606683729504
```

Final density-matrix difference versus v0.16:

```text
7.09620610556202e-15
```

The v0.17 result therefore reproduces the accepted v0.16 trajectory to floating-point
precision after the online controller relaxes the graph.

## Snapshot convergence

Edge-score convergence:

```text
S monotone:
True

H monotone:
True

finest S error:
0.00026621461118036714

finest H error:
0.00015704048353805722
```

Local-score-budget convergence:

```text
S monotone:
True

H monotone:
True

strictest S error:
0.0

strictest H error:
0.0
```

## Bounded-locality scaling

At `N=160`:

```text
active edges:
317

all possible off-diagonal edges:
12720

exact S/H/T score checks:
630

pair factorizations:
790

dense canonical pairs:
12880

pair reduction:
93.866 %

dense matrix assembly:
3.823627 s

sparse matrix assembly:
0.029290 s

diagnostic assembly speedup:
130.544 x
```

Fitted exponents over `N=20,40,80,160`:

```text
active edges:
1.0324491201728092

KD-tree spatial candidates:
1.0557080719105292

exact S/H/T checks:
1.0557080719105292

pair factorizations:
1.043904194562303

dense canonical pairs:
1.97980963642741
```

These exponents apply only to the deterministic bounded-locality chain.

## Complexity ledger

```text
dense audits:
8

audit pair factorizations:
506

audit time:
0.148051 s

score relaxations:
1

search-floor relaxations:
1

propagation pair factorizations:
15730

candidate pair factorizations:
88

graph/S-H-T score time:
8.550292 s

sparse matrix assembly:
0.495324 s

sparse T assembly:
1.460446 s

sparse Cayley solves:
0.085747 s

TDSE-defect work:
0.481460 s

total adaptive run:
13.552188 s
```

The periodic dense audit intentionally reintroduces an $O(N^2)$ correctness cost.
v0.17 therefore does not claim end-to-end asymptotically sparse dynamics.

## Representative examples

The following examples were executed successfully during validation:

```text
examples/71_v017_edge_importance.py
examples/72_v017_online_audit.py
examples/73_v017_threshold_convergence.py
examples/74_v017_budget_convergence.py
examples/75_v017_scaling.py
examples/76_v016_v017_comparison.py
```

`examples/77_recompute_v017_campaign.py` is the complete campaign recomputation wrapper.

## PySCF status

The inherited PySCF/state-tracking/gauge-graph infrastructure remains in the cumulative
repository.

v0.17 does not claim a calibrated molecular S/H/T edge-error model or a production
PySCF sparse residual controller.

See `V17_PYSCF_ERROR_CONTROL_BRIDGE.md`.
