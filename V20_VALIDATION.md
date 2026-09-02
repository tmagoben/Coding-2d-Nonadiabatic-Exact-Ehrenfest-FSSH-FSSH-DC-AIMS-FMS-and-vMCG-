# v0.20 Validation Contract

v0.20 is accepted only if the molecular sparse machinery passes moving-dynamics,
matrix-convergence, scaling, indexed-cache, and online-controller checks.

## Canonical sparse dynamics

Required:

```text
dense-metric coefficient error <= 1e-3
center error <= 1e-12
norm drift <= 1e-10
average matrix sparsity >= 75%
new electronic-point reduction >= 40%
sampled audit failures = 0
```

Measured:

```text
coefficient error: 0.000670432070175427
center error: 0.0
norm drift: 8.881784197001252e-16
average sparsity: 81.05 %
electronic-point reduction: 49.30 %
sampled audit failures: 0
```

All pass.

## Final dense sentinel

Limits:

```text
S error <= 0.005
H error <= 0.003
T seed error <= 0.020
```

Measured:

```text
S: 0.002954918852389199
H: 0.0016775567668417452
T: 0.01718013850955588
```

All pass.

## Score-threshold convergence

| Enter score | Active edges | S error | H error | T error |
|---:|---:|---:|---:|---:|
| 0.120 | 21 | 0.019303621 | 0.0029620109 | 0.052087976 |
| 0.080 | 21 | 0.019303621 | 0.0029620109 | 0.052087976 |
| 0.050 | 21 | 0.019303621 | 0.0029620109 | 0.052087976 |
| 0.030 | 23 | 0.011138901 | 0.0025357232 | 0.042026352 |
| 0.020 | 23 | 0.011138901 | 0.0025357232 | 0.042026352 |
| 0.010 | 27 | 0.0070879361 | 0.0020982307 | 0.029257706 |
| 0.005 | 34 | 0.0034262757 | 0.001697729 | 0.017544481 |

The S/H/T errors are all nonincreasing. The finest row also satisfies the release
matrix-error limits.

## Local-score-budget convergence

| Budget | Active edges | Promoted | Score L2 | S error | H error | T error |
|---:|---:|---:|---:|---:|---:|---:|
| 1e+09 | 23 | 0 | 0.037174371 | 0.011138901 | 0.0025357232 | 0.042026352 |
| 0.05 | 23 | 0 | 0.037174371 | 0.011138901 | 0.0025357232 | 0.042026352 |
| 0.02 | 29 | 6 | 0.019800169 | 0.005908745 | 0.0020486794 | 0.028194466 |
| 0.01 | 36 | 13 | 0.0098800342 | 0.002953481 | 0.0016790767 | 0.017194256 |
| 0.005 | 53 | 30 | 0.0048760583 | 0.0014519635 | 0.00071056245 | 0.010443794 |
| 0 | 85 | 62 | 0 | 6.1534034e-07 | 4.5285153e-07 | 6.9841207e-06 |

At zero budget:

```text
S = 6.153403390687608e-07
H = 4.528515298636401e-07
T = 6.984120676640802e-06
```

all lie below the strict zero-budget acceptance limits.

## Bounded-locality scaling

Required exponents:

```text
active <= 1.20
pair checks <= 1.10
provider misses <= 1.10
dense pairs >= 1.90
```

Measured:

```text
active = 1.117285787300988
pair checks = 1.0324491201728092
provider misses = 1.0213244854837158
dense pairs = 2.0213244854837162
```

At $N=160$ the pair-check reduction is

```text
97.51 %
```

and exceeds the required 95%.

## Online search controller

The deliberately aggressive test begins at search floor `0.90`.

First audit:

```text
maximum sampled score = 0.45067614545809304
passed = False
```

The controller changes the search floor to

```text
0.45
```

and the immediate re-audit passes.

## Indexed cache

Automated tests verify exact nearest-neighbor agreement against brute force for random
multidimensional points.

The molecular provider regression also verifies the indexed provider preserves tracked
electronic behavior.

## PySCF status

```text
installed: False
runtime validated: False
```

No real PySCF sparse trajectory is claimed.

## Release result

```text
passed = True
```
