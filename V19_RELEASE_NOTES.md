# v0.19 Release Notes

v0.19 is the **molecular/direct-dynamics integration release**.

## New modules

```text
molecular_snapshot_v19.py
analytic_molecular_backend_v19.py
state_tracking_v19.py
molecular_direct_provider_v19.py
molecular_gauge_graph_v19.py
molecular_direct_dynamics_v19.py
pyscf_molecular_bridge_v19.py
v19_benchmark.py
```

## Main additions

```text
overlap-capable molecular electronic snapshot contract
Cartesian -> generalized electronic projection validation
nearest-anchor tracking for branched center/centroid requests
exact geometry cache
provider cost estimates
explicit backend retry/fallback policy
polynomial maximum-overlap state assignment
center-centroid electronic gauge graph
provider-neutral spawned molecular Gaussian dynamics
raw PySCF SA-CASSCF snapshot adapter
many-electron overlap-engine bridge
```

## Major algorithmic change

State assignment changes from

$$
O(n_s!)
$$

permutation enumeration to

$$
O(n_s^3)
$$

for the best Hungarian assignment and

$$
\boxed{O(n_s^4)}
$$

when retaining the exact second-best ambiguity margin.

The release validates this path for

```text
nstate = 16
```

while small-state tests verify exact agreement with the older exhaustive method.

## Molecular projection/tracking

```text
maximum tracked energy error:
5.204170427930421e-18

maximum tracked gradient error:
5.204170427930421e-18

maximum tracked NAC error:
4.440892098500626e-16

untracked scrambled-root energy error:
0.03542632418556527
```

## Molecular gauge graph

```text
S clean/scrambled difference:
0.0

H clean/scrambled difference:
0.0

H Hermiticity error:
1.2668532396679003e-18
```

## Spawned direct dynamics

```text
spawn events:
1

coefficient difference vs analytic reference:
3.3009808712653696e-16

center difference:
4.47545209131181e-16

maximum norm drift:
7.771561172376096e-16
```

## PySCF status

```text
installed in build environment:
False

runtime validated:
False
```

No real PySCF calculation is claimed.

## Release status

```text
passed = True
```

The cumulative automated regression suite reports:

```text
224 passed
```

Full details are recorded in `V19_BUILD_VALIDATION.md`.

## Scope

v0.19 is still a research prototype. It is not production AIMS, does not yet use a
sparse molecular pair-centroid graph for large bases, and does not yet include SOC.
