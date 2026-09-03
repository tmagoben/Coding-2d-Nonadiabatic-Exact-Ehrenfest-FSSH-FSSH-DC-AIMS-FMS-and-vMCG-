# v0.19 Build Validation Report

Validated on 2026-08-13.

## Source validation

```text
313 Python files parsed successfully with Python AST.
```

## Cumulative automated regression suite

```text
224 passed in 14.09 s
```

The test suite remains cumulative from v0.1 through v0.19.

New v0.19 coverage includes:

```text
Cartesian molecular -> generalized coordinate projection
root permutation/sign-scramble recovery
nearest-anchor state tracking
nonsequential query order after reference seeding
exact geometry cache reuse
provider cost classification
bounded explicit backend fallback
PySCF-like injected overlap-engine path
center-centroid molecular gauge graph
gauge-invariant S/H matrix construction
short spawned molecular direct dynamics
polynomial state assignment
exact small-manifold agreement with exhaustive assignment
16-state tracking diagnostic
complete v0.19 release campaign acceptance
```

## Canonical campaign

Machine-readable result:

```text
results/v019_molecular_direct_dynamics_campaign.json
```

Release status:

```text
passed = True
```

All configured acceptance checks pass.

## Molecular Cartesian projection

The deterministic validation backend embeds the exact 2D LVC problem into a two-atom
Cartesian geometry.

Measured maximum tracked errors:

```text
energy:
5.204170427930421e-18

gradient:
5.204170427930421e-18

NAC:
4.440892098500626e-16
```

The deliberately untracked scrambled-root backend reaches:

```text
raw energy-order error:
0.03542632418556527
```

so the state-tracking validation is nontrivial.

## Nonsequential tracking

After one reference geometry seeds the electronic labels:

```text
shuffled-order max energy error:
5.204170427930421e-18

shuffled-order max NAC error:
4.440892098500626e-16

tracking ambiguities:
0
```

The branched-query stress test passes.

## Molecular gauge graph

```text
basis size:
3

nodes:
6

edges:
6

clean/scrambled S difference:
0.0

clean/scrambled H difference:
0.0

S Hermiticity error:
0.0

H Hermiticity error:
1.2668532396679003e-18

condition number:
7.385474223582431
```

## Spawned direct dynamics

All three runs

```text
analytic generalized provider
clean molecular Cartesian provider
scrambled-root molecular provider + tracking
```

produce the same one-child spawn event.

Measured analytic-vs-scrambled differences:

```text
coefficient difference:
3.3009808712653696e-16

center difference:
4.47545209131181e-16

momentum difference:
6.646852490838693e-13

generalized norm drift:
7.771561172376096e-16
```

Final center-centroid graph audit:

```text
{'H_hermiticity_error': 0.0, 'S_hermiticity_error': 0.0, 'basis_size': 2, 'condition_number': 1.0001619386564635, 'graph_edges': 2, 'graph_nodes': 3}
```

## Cache and failure behavior

Scrambled scan provider:

```text
cache hits:
17

cache misses:
23

backend attempts:
23

tracking ambiguities:
0
```

Failure test:

```text
backend failures:
1

fallback uses:
1

fallback distance:
0.039999999999999994
```

The default provider failure policy remains `raise`. The validated fallback is explicit
and opt-in.

## Provider-cost interface

```text
cached normalized cost:
0.1

nearby normalized cost:
0.5

new normalized cost:
5.0
```

These values are deterministic scheduling units, not measured PySCF costs.

## State-assignment scaling

The older exhaustive assignment has factorial growth.

v0.19 validates:

```text
best assignment:
O(nstate^3)

best + exact second-best ambiguity margin:
O(nstate^4)
```

Large-state diagnostic:

```text
nstate:
16

valid permutation:
True

best score:
3.9437377812776653

second-best score:
3.9410855307265416

diagnostic wall time:
0.00030213800005185476 s
```

Automated tests separately verify exact agreement with the exhaustive implementation for
2, 3, 4, and 5 states.

## Representative examples executed

```text
examples/87_v019_molecular_projection.py
examples/88_v019_tracking_order.py
examples/89_v019_molecular_gauge_graph.py
examples/90_v019_direct_dynamics.py
examples/91_v019_cache_failure_cost.py
examples/92_v019_state_assignment.py
examples/93_v019_pyscf_bridge.py
```

## PySCF status

```text
PySCF installed in build environment:
False

real PySCF runtime validated:
False
```

The following bridge is present and regression-tested at the interface level:

```text
PySCFRawSnapshotBackendV19
        +
pyscf_snapshot_overlap_engine_v19
        +
TrackedMolecularDirectProviderV19
```

No real PySCF trajectory or ab-initio numerical result is claimed.

## Scientific limitations

v0.19 remains a research prototype.

The current release does not yet provide:

```text
production AIMS matrix elements
large-basis sparse molecular centroid construction
dynamic spatial indexing of a very large electronic cache
asynchronous electronic-structure scheduling
real PySCF runtime validation in this environment
complex electronic/NAC data contracts
spin-orbit coupling dynamics
```

These limitations are explicit in the theory and PySCF protocol documents.
