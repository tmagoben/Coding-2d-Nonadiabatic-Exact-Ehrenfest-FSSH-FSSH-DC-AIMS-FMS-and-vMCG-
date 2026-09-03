# v0.20 Benchmark Results

Canonical output:

```text
results/v020_sparse_molecular_campaign.json
```

## Moving sparse-vs-dense dynamics

```text
possible pairs: 190
active edges: 36
average sparsity: 81.05 %

dense-metric coefficient error: 0.000670432070175427
center error: 0.0
momentum error: 0.0
norm drift: 8.881784197001252e-16
```

## Final dense sentinel

```text
S error: 0.002954918852389199
H error: 0.0016775567668417452
T seed error: 0.01718013850955588
```

## Electronic work

```text
sparse production cache misses: 2149
dense reference cache misses: 4239
reduction: 49.30 %
```

Diagnostic wall time:

```text
sparse: 3.945 s
dense: 7.105 s
speedup: 1.80 x
```

## Molecular threshold convergence

| Enter | Active | S error | H error | T error |
|---:|---:|---:|---:|---:|
| 0.120 | 21 | 0.019303621 | 0.0029620109 | 0.052087976 |
| 0.080 | 21 | 0.019303621 | 0.0029620109 | 0.052087976 |
| 0.050 | 21 | 0.019303621 | 0.0029620109 | 0.052087976 |
| 0.030 | 23 | 0.011138901 | 0.0025357232 | 0.042026352 |
| 0.020 | 23 | 0.011138901 | 0.0025357232 | 0.042026352 |
| 0.010 | 27 | 0.0070879361 | 0.0020982307 | 0.029257706 |
| 0.005 | 34 | 0.0034262757 | 0.001697729 | 0.017544481 |

## Molecular local-budget convergence

| Budget | Active | Promoted | Score L2 | S error | H error | T error |
|---:|---:|---:|---:|---:|---:|---:|
| 1e+09 | 23 | 0 | 0.037174371 | 0.011138901 | 0.0025357232 | 0.042026352 |
| 0.05 | 23 | 0 | 0.037174371 | 0.011138901 | 0.0025357232 | 0.042026352 |
| 0.02 | 29 | 6 | 0.019800169 | 0.005908745 | 0.0020486794 | 0.028194466 |
| 0.01 | 36 | 13 | 0.0098800342 | 0.002953481 | 0.0016790767 | 0.017194256 |
| 0.005 | 53 | 30 | 0.0048760583 | 0.0014519635 | 0.00071056245 | 0.010443794 |
| 0 | 85 | 62 | 0 | 6.1534034e-07 | 4.5285153e-07 | 6.9841207e-06 |

## Bounded-locality scaling

| N | Active | Pair checks | Dense pairs | Reduction | Electronic points |
|---:|---:|---:|---:|---:|---:|
| 20 | 15 | 37 | 190 | 80.53% | 57 |
| 40 | 32 | 77 | 780 | 90.13% | 117 |
| 80 | 71 | 157 | 3160 | 95.03% | 237 |
| 160 | 152 | 317 | 12720 | 97.51% | 477 |

Fitted exponents:

```text
active = 1.117285787300988
pair checks = 1.0324491201728092
electronic points = 1.0213244854837158
dense pairs = 2.0213244854837162
```

At `N=160`:

```text
active edges = 152
pair checks = 317
dense pairs = 12720
pair-check reduction = 97.51 %
matrix sparsity = 98.81 %
```

## Online controller

```text
initial missed-edge score = 0.45067614545809304
search floor: 0.90 -> 0.45
final audit passed = True
```

## Release status

```text
passed = True
```
