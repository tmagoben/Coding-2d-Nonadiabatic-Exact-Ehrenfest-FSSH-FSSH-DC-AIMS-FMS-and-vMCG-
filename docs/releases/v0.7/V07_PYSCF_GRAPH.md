# v0.7 PySCF Gauge-Graph Workflow

v0.7 reuses the explicit SA-CASSCF snapshot and many-electron overlap machinery from
v0.6.

## 1. Node data

Each electronic graph node can be a

```text
TBF center
pair centroid
additional overlap-check point
```

with a `CASSCFWavefunctionSnapshot` containing

```text
PySCF Mole
MO coefficients
CI roots
ncore
ncas
nelecas
```

## 2. Edge construction

For selected neighboring nodes $u,v$,

```python
O_uv = casscf_state_overlap_matrix(snapshot_u, snapshot_v)
```

is evaluated with the full nonorthogonal core+active overlap engine from v0.6.

Then

```python
graph.add_overlap(u, v, O_uv)
```

projects the overlap to its nearest unitary link.

## 3. Convenience builder

```python
from gaussian_dynamics.pyscf_gauge_graph import build_snapshot_gauge_graph

graph = build_snapshot_gauge_graph(
    snapshots,
    edge_pairs,
)
```

where `snapshots` is a mapping

```text
node_id -> CASSCFWavefunctionSnapshot
```

and `edge_pairs` declares which electronic calculations should be compared.

## 4. TBF / centroid connectivity

For TBF nodes

```text
t0, t1, t2
```

and centroid nodes

```text
c01, c02, c12
```

the helper

```python
tbf_centroid_edge_pairs(...)
```

returns

```text
t0--c01--t1
t0--c02--t2
t1--c12--t2
```

which contains a closed gauge loop.

## 5. Diagnostics to inspect before propagation

For every edge, inspect

```text
singular values of O_uv
||O_uv^dagger O_uv - I||_F
```

Large singular-value loss or a large unitarity defect suggests that

- the geometry separation is too large;
- the selected state manifold is incomplete;
- roots have changed character strongly;
- the graph needs intermediate nodes.

For every fundamental cycle, inspect

```python
W = graph.wilson_loop(cycle)
```

and record

```text
eigenvalues(W)
trace(W)
```

A nontrivial Wilson matrix may represent genuine geometric holonomy and should not be
forced to the identity simply because a smoother gauge is desired.

## 6. Relationship to v0.6 sequential tracking

v0.6 remains appropriate for one ordered trajectory path.

v0.7 is intended when the electronic calculations form a network.

The two layers are complementary:

```text
v0.6
sequential root identity + many-electron overlap

v0.7
network connection + loops + global gauge consistency
```

## 7. Recommended ab initio workflow

```text
1. Run SA-CASSCF at required TBF centers and pair centroids.
2. Save every CASSCFWavefunctionSnapshot.
3. Connect nearby nodes by many-electron overlap edges.
4. Check edge singular values/unitarity defects.
5. Compute fundamental Wilson loops.
6. Build a tree gauge for transparent diagnostics.
7. Optionally synchronize the graph to distribute finite-step noise.
8. Transport TBF electronic coefficient vectors to pair-centroid gauges.
9. Construct graph-covariant pair electronic factors.
10. Only then assemble Gaussian pair matrices.
```

## 8. Important limitation

The current graph builder assumes that every node contains the same number of tracked
roots/subspace dimension.

A production code that dynamically changes the electronic state manifold would need a
rectangular-subspace matching layer and explicit rules for enlarging/reducing the
manifold.
