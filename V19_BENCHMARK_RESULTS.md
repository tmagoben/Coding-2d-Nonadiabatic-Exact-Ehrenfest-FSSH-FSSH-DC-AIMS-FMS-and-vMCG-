# v0.19 Benchmark Results

Canonical output:

```text
results/v019_molecular_direct_dynamics_campaign.json
```

## Electronic data projection/tracking

```text
scan points:
17

maximum clean/tracked energy error:
5.204170427930421e-18

maximum clean/tracked gradient error:
5.204170427930421e-18

maximum clean/tracked NAC error:
4.440892098500626e-16

raw untracked scrambled energy error:
0.03542632418556527
```

The large raw error confirms that state tracking is essential.

## Nonsequential query test

```text
shuffled maximum energy error:
5.204170427930421e-18

shuffled maximum NAC error:
4.440892098500626e-16

ambiguities:
0
```

## Molecular gauge graph

```text
basis size:
3

nodes:
6

edges:
6

S clean/scrambled difference:
0.0

H clean/scrambled difference:
0.0

H Hermiticity error:
1.2668532396679003e-18

condition number:
7.385474223582431
```

## Direct dynamics

Spawn event:

```text
step:
1

time:
0.02

target state:
0
```

Analytic-vs-scrambled molecular differences:

```text
coefficient difference:
3.3009808712653696e-16

center difference:
4.47545209131181e-16

momentum difference:
6.646852490838693e-13

maximum norm drift:
7.771561172376096e-16
```

## Cache/failure contract

Provider scan/cache:

```text
cache hits:
17

cache misses:
23

backend attempts:
23
```

Failure demonstration:

```text
backend failures:
1

fallback uses:
1

fallback distance:
0.039999999999999994
```

## Cost estimate

```text
cached:
0.1

nearby:
0.5

new:
5.0
```

These are normalized scheduling costs, not ab-initio timings.

## State assignment

```text
nstate:
16

valid permutation:
True

best score:
3.9437377812776653

second-best score:
3.9410855307265416

diagnostic time:
0.00030213800005185476 s
```

Algorithmic complexity:

```text
legacy: O(n!)
v0.19: O(n^4) with exact second-best margin
```

## PySCF

```text
installed in build environment:
False

runtime validated:
False

bridge:
PySCFRawSnapshotBackendV19 + casscf_state_overlap_matrix
```

## Acceptance

```text
passed = True
```

All configured checks pass.
