# v0.19 Molecular Direct-Dynamics Benchmark Campaign

## Purpose

The v0.19 campaign asks whether the existing Gaussian dynamics machinery can receive
molecular-style Cartesian electronic-structure data without losing state identity,
gauge consistency, or reproducibility.

The campaign is intentionally deterministic and independently checkable.

## Synthetic molecular model

A two-atom H2-like Cartesian geometry is used with two orthonormal collective modes:

```text
mode 0: relative x stretch-like displacement
mode 1: relative transverse y displacement
```

The electronic Hamiltonian is the exact 2-state analytic LVC conical-intersection model.

This gives an exact generalized-coordinate reference while forcing the software through

```text
MolecularGeometry
CartesianElectronicStructurePoint
Cartesian gradients
Cartesian NACs
atomic masses
LinearGeometryMap
generalized-coordinate projection
```

## Root-scrambling stress test

The raw backend deliberately changes:

```text
state order
state signs
```

with geometry.

The tracked provider must recover the same generalized electronic data as the clean
analytic reference.

## Scan

```text
17 generalized-coordinate points
x from -0.8 to 0.8
y near 0.35
```

The first point establishes the tracked reference labels.

A second provider then evaluates the remaining points in nonsequential order after the
same reference seed.

## Gauge graph

A 3-Gaussian basis is used.

The graph contains:

```text
3 center nodes
3 pair-centroid nodes
6 center-centroid edges
```

Static gauge-covariant S/H matrices are compared between clean and root-scrambled
molecular backends.

## Direct dynamics

Initial condition:

```text
one upper-surface Gaussian
q = (-0.45, 0.35)
p = (20, 0)
width = 0.8 I
```

Propagation:

```text
dt = 0.02
steps = 50
spawn threshold = 0.001
maximum basis = 2
```

The energy-conserving child is created at the first step.

Three runs are compared:

```text
analytic generalized provider
clean molecular provider
scrambled + tracked molecular provider
```

## Failure/cost contracts

The campaign also validates:

```text
exact cache hits
near-cache cost estimate
new-point cost estimate
explicit backend failure
bounded nearest-cache fallback
```

## State-tracking scaling

A deterministic 16-state overlap assignment validates the Hungarian v0.19 path.

Small-state automated tests compare the exact best and second-best scores directly with
the older exhaustive implementation.

## PySCF

PySCF is not required for campaign execution.

The real backend bridge is present but runtime validation is deferred to an environment
with PySCF installed.

## Reproduction

```bash
python examples/94_recompute_v019_campaign.py
```

Canonical output:

```text
results/v019_molecular_direct_dynamics_campaign.json
```
