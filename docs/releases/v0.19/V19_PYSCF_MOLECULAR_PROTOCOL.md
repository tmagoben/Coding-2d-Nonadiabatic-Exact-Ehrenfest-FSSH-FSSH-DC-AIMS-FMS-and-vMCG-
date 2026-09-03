# v0.19 PySCF Molecular Direct-Dynamics Protocol

This document describes how the v0.19 architecture is intended to be used with the
existing PySCF SA-CASSCF machinery.

It is a protocol, not a runtime-validation report.

PySCF was not installed in the v0.19 build environment.

## 1. Use the raw snapshot adapter

Do not use the sequential tracked backend as the electronic provider for a branched
Gaussian basis.

Use:

```python
PySCFRawSnapshotBackendV19(config)
```

The adapter exposes one raw SA-CASSCF point plus one
`CASSCFWavefunctionSnapshot`.

State tracking is then owned by:

```python
TrackedMolecularDirectProviderV19
```

with

```python
overlap_engine=pyscf_snapshot_overlap_engine_v19
```

This separates the electronic calculation from branched state/gauge tracking.

## 2. Why sequential tracking is insufficient

A sequential geometry path has a natural previous point. A Gaussian basis generally
does not.

Electronic requests may arrive as TBF centers and pair centroids in an order that is
set by software, not physical time. v0.19 instead tracks a new point against the nearest
trusted cached electronic anchor.

## 3. Establish the reference frame first

The first accepted geometry defines tracked labels and phases. Therefore a molecular
run should explicitly evaluate the physical initial geometry before any branched
centroid requests.

Recommended ordering:

```text
initial molecular geometry
        ↓
accept reference electronic snapshot
        ↓
initial TBF centers
        ↓
pair centroids / spawned centers
```

## 4. PySCF NAC convention

The internal dynamics convention remains

$$
d_{ij}
=
\langle\Phi_i|\nabla_R\Phi_j\rangle.
$$

The inherited PySCF backend interprets the state tuple as:

```text
state = (ket, bra)
```

Therefore internal `d[i,j]` is requested with:

```text
state = (j, i)
mult_ediff = False
```

The `mult_ediff=True` quantity remains diagnostic only.

## 5. Many-electron cross-geometry overlap

The v0.6 overlap engine computes

$$
O_{IJ}
=
\langle
\Psi_I(R_A)
|
\Psi_J(R_B)
\rangle.
$$

Its ingredients are cross-AO overlap, core+active MO cross overlap, embedded
core+active CI vectors, and the PySCF FCI overlap with the nonorthogonal orbital metric.

The resulting overlap matrix is passed to the v0.19 polynomial state assignment.

## 6. State assignment

v0.19 no longer enumerates all root permutations.

For `nstate` roots:

```text
best assignment:
O(nstate^3)

best + exact second-best ambiguity margin:
O(nstate^4)
```

This is still intended for modest state manifolds, but is substantially more scalable
than factorial enumeration.

## 7. Ambiguity policy

Recommended production-like validation settings are:

```text
minimum assigned overlap: system dependent
minimum score margin: system dependent
ambiguity policy: raise
```

Near real degeneracies, root-by-root assignment may become intrinsically unstable. The
existing subspace/Procrustes infrastructure should then be used rather than forcing a
unique individual-state identity.

## 8. Geometry cache

Exact cached generalized coordinates return without rerunning PySCF.

Warm-start orbitals may reduce SCF/CASSCF initialization cost, but they do not solve
state identity, phase/gauge identity, or many-electron overlap continuity.

## 9. Cost-model calibration

The current normalized costs are placeholders.

A real PySCF calibration should record at least:

```text
SCF wall time
SA-CASSCF wall time
gradient wall time
NAC wall time
cross-geometry overlap wall time
state tracking wall time
```

## 10. Failure policy

The default should remain:

```text
failure_policy = "raise"
```

The opt-in nearest-cache fallback is intended only for controlled robustness
experiments. It should use a small explicit distance bound, mark stale reuse in
metadata, and never promote the fallback point to a trusted gauge anchor.

## 11. Sparse molecular graph

The v0.19 validation graph builds all pair centroids. A large molecular calculation
should not.

Instead, centroid electronic calculations should be requested only for active nuclear
Gaussian locality edges:

$$
P\ll N^2.
$$

The intended future architecture is:

```text
nuclear sparse locality graph
        ↓
active Gaussian pair list
        ↓
only required molecular pair-centroid electronic points
        ↓
electronic gauge links
        ↓
sparse molecular S/H/T construction
```

## 12. Validation sequence for a real molecule

Before accepting a molecular direct-dynamics calculation, validate:

```text
1. SCF convergence
2. SA-CASSCF convergence
3. state energies
4. gradients
5. NAC convention/sign sanity
6. cross-geometry state overlaps
7. tracking assignment confidence
8. gauge-loop diagnostics where relevant
9. generalized-coordinate mass matrix
10. Gaussian overlap conditioning
11. timestep convergence
12. basis/spawn convergence
13. electronic-cache sensitivity
```

Only after these are stable should SOC be added.

## 13. Current release boundary

v0.19 provides the architecture and deterministic regression tests for this workflow.

It does not claim that a real PySCF trajectory was executed in the build environment.
