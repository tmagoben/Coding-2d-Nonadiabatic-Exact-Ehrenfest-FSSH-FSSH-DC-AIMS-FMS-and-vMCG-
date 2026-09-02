# v0.8 PySCF Dynamic-Graph Workflow

The PySCF-facing v0.8 workflow deliberately avoids forcing all electronic evaluations through one sequential state tracker.

## 1. Raw point plus CASSCF snapshot

Use

```python
point, snapshot = backend.evaluate_raw_with_snapshot(geometry)
```

from `PySCFTrackedSACASSCFBackend`.

Despite the class name inherited from v0.6, this method performs one raw SA-CASSCF evaluation and returns its wavefunction snapshot without applying the sequential root-tracking history.

That is the preferred object for graph dynamics.

## 2. Incremental graph

Create

```python
builder = IncrementalSnapshotGaugeGraph(nstates)
```

and add the first point:

```python
builder.add_cartesian_point(
    node_id,
    snapshot,
    point,
)
```

For a new TBF-center or centroid node, specify only the existing nodes to which a many-electron overlap edge is needed:

```python
builder.add_cartesian_point(
    new_node,
    new_snapshot,
    new_point,
    connect_to=[parent_node, centroid_neighbor],
)
```

The builder computes only those overlap matrices.

## 3. Why raw snapshots are preferable on a graph

A sequential root tracker assumes

```text
R0 -> R1 -> R2 -> ...
```

but a spawned Gaussian graph can request

```text
parent center
child center
pair centroid
other pair centroid
```

in an order that is not one physical path.

The graph should therefore determine relative electronic gauge using explicit overlap edges rather than relying on mutable call order.

## 4. Cartesian mode

`add_cartesian_point` uses all $3N$ Cartesian coordinates directly:

```text
gradients_cart : (nstate,natom,3)
nac_cart       : (nstate,nstate,natom,3)
```

which are flattened to the derivative-Hamiltonian field consumed by the graph registry.

If a reduced coordinate system is desired, project first with the v0.5 `LinearGeometryMap` machinery and add generalized operator matrices instead.

## 5. Raw-overlap quality

The graph edge ultimately stores a unitary polar link, but always inspect the raw overlap diagnostics:

```python
builder.edge_diagnostics()
```

including

$$
\sigma_k(O)
$$

and

$$
\|O^\dagger O-I\|_F.
$$

A perfect-looking unitary polar factor does not mean the original selected electronic subspaces had strong overlap.

## 6. PySCF runtime status

PySCF is optional and is not installed in the release build environment. The real backend path therefore remains source/API aligned and fake-backend regression tested, while all dynamic graph mathematics is runtime validated on the analytic CI provider.
