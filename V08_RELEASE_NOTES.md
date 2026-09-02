# v0.8 Release Notes

## Theme

**Time-dependent graph-AIMS-style propagation.**

## New modules

```text
temporal_electronic.py
dynamic_gauge_graph.py
moving_graph_gaussian.py
dynamic_graph_aims.py
incremental_snapshot_graph.py
```

## Main scientific additions

- temporal overlap/local-diabatic electronic propagation;
- direct comparison to explicit-NAC propagation;
- incremental time-labelled gauge graph;
- TBF-center and pair-centroid nodes created during propagation;
- metric-compatible moving-basis $T$ matrix;
- dynamic zero-amplitude spawning;
- coupled child amplitude growth;
- public raw PySCF point+snapshot entry point;
- incremental many-electron PySCF graph construction.
