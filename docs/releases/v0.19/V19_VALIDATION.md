# v0.19 Validation Contract

v0.19 separates molecular integration validation into six layers:

1. Cartesian/generalized coordinate correctness;
2. electronic root/gauge tracking;
3. molecular gauge-graph covariance;
4. direct-dynamics equivalence;
5. caching/failure behavior;
6. state-assignment scaling.

## 1. Cartesian projection

Thresholds:

```text
maximum energy error     <= 1e-11
maximum gradient error   <= 1e-11
maximum NAC error        <= 1e-11
```

Measured:

```text
energy:
5.204170427930421e-18

gradient:
5.204170427930421e-18

NAC:
4.440892098500626e-16
```

All pass.

## 2. Tracking must correct a real failure mode

The deliberately scrambled raw backend has maximum untracked energy-order error

```text
0.03542632418556527
```

which is required to exceed `1e-4`.

Thus the tracking test is not vacuous.

After tracking the error falls to machine precision.

## 3. Nonsequential request order

After establishing one reference seed, geometries are requested in a deliberately
shuffled order.

Measured:

```text
maximum energy error:
5.204170427930421e-18

maximum NAC error:
4.440892098500626e-16

tracking ambiguities:
0
```

All pass.

This validates order tolerance after a reference electronic frame has been established.

## 4. Gauge-graph matrices

Thresholds:

```text
clean/scrambled S difference <= 1e-10
clean/scrambled H difference <= 1e-10
S/H Hermiticity error        <= 1e-10
```

Measured:

```text
S difference:
0.0

H difference:
0.0

S Hermiticity:
0.0

H Hermiticity:
1.2668532396679003e-18

condition number:
7.385474223582431
```

All pass.

## 5. Spawned direct dynamics

The benchmark compares:

```text
analytic generalized LVC dynamics
clean molecular Cartesian bridge
scrambled molecular Cartesian bridge + tracking
```

Each produces exactly one spawn event.

Measured analytic-vs-scrambled differences:

```text
coefficients:
3.3009808712653696e-16

centers:
4.47545209131181e-16

momenta:
6.646852490838693e-13

generalized norm drift:
7.771561172376096e-16
```

The coefficient, center, and norm acceptance criteria all pass.

## 6. Cache behavior

The provider scan records:

```text
cache hits:
17

cache misses:
23

tracking ambiguities:
0
```

Required:

```text
cache hits >= 10
tracking ambiguities = 0
```

Both pass.

## 7. Failure policy

The failure benchmark performs one successful electronic calculation and then requests
a point inside a deterministic backend failure region.

The failed point is within the configured fallback radius.

Measured:

```text
backend failures:
1

fallback uses:
1

fallback distance:
0.039999999999999994
```

Exactly one fallback occurs.

The default policy remains `raise`; this fallback is opt-in.

## 8. State assignment

v0.19 validates the polynomial assignment algorithm at

```text
nstate = 16
```

with a valid permutation.

The best/second-best assignment scores are both finite and the full small-manifold test
suite verifies agreement with the old exhaustive algorithm for 2-5 states.

## 9. PySCF scope

```text
PySCF installed:
False

real PySCF runtime validated:
False
```

The snapshot bridge and injected overlap-engine path are tested without requiring a real
PySCF calculation.

No real PySCF numerical result is claimed.

## 10. Release result

```text
passed = True
```

All configured v0.19 checks pass.
