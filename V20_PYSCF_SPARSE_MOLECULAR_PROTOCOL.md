# v0.20 PySCF Sparse Molecular Protocol

v0.20 completes the sparse molecular architecture, but PySCF is not installed in the
build environment.

## Provider stack

Use:

```text
PySCFRawSnapshotBackendV19
        ↓
pyscf_snapshot_overlap_engine_v19
        ↓
IndexedTrackedMolecularDirectProviderV20
```

The v0.20 provider adds indexed nearest-anchor tracking to the v0.19 many-electron
snapshot bridge.

## Do not construct every pair centroid

For $N$ TBFs, do not request all $N(N-1)/2$ molecular centroids by default.

Use:

```text
Gaussian position/width KD-tree
        ↓
pair-specific nuclear-overlap screen
        ↓
only admitted candidate centroids
        ↓
molecular S/H/T score
        ↓
active sparse edges
```

For bounded locality, candidate electronic calculations should scale with local graph
degree rather than $N^2$.

## Candidate centroid cost is real electronic cost

Every admitted molecular candidate may require SCF, SA-CASSCF, gradients, NACs, and a
wavefunction snapshot unless already cached.

Therefore `search_overlap_floor` is an electronic-cost control and must be calibrated
with sampled omitted-edge audits.

## Sampled audits

Recommended samples combine:

```text
near-boundary priority pairs
random omitted pairs
known difficult electronic regions
recent failure/retry neighborhoods
```

If a sampled omitted edge has exact score above the enter threshold, lower the
geometric search floor and rebuild.

Do not relax physical accuracy thresholds merely because the electronic calculation is
expensive.

## Dense validation must use separate caches

If a complete dense validation snapshot is feasible, run it through an independent
provider/cache.

Otherwise the audit precomputes every centroid and invalidates sparse electronic-cost
measurements.

For realistic systems, full dense sentinels should usually be restricted to small
validation molecules or short selected trajectory windows.

## Indexed tracking cache

Exact geometry hits remain hash lookups.

New geometries use a buffered KD-tree nearest trusted anchor plus many-electron
cross-geometry state overlap.

The physical initial geometry should seed the electronic label/gauge frame before
branched center/centroid queries begin.

## PySCF NAC convention

The internal convention remains

$$
d_{ij}=\langle\Phi_i|\nabla_R\Phi_j\rangle.
$$

PySCF's state tuple is interpreted as `(ket, bra)`, so internal `d[i,j]` uses

```text
state = (j, i)
mult_ediff = False
```

`mult_ediff=True` remains diagnostic only.

## Sparse linear algebra

The current v0.20 molecular matrices are scalar TBF matrices in the discrete
pair-centroid approximation.

A future spinor/SOC extension should preserve the same active-edge sparse graph and
replace scalar entries with electronic blocks rather than reverting to dense
all-pairs construction.

## Real molecular validation before SOC

A real PySCF validation should demonstrate:

```text
state tracking continuity
new electronic-point reduction
matrix sentinel agreement on small snapshots
timestep convergence
basis convergence
sparse-threshold convergence
sampled-audit stability
cache/failure reproducibility
```

Only then should SOC be layered onto the molecular sparse engine.

## Current claim boundary

v0.20 validates the sparse machinery using the deterministic Cartesian LVC molecular
backend. No real PySCF sparse trajectory is claimed in this environment.
