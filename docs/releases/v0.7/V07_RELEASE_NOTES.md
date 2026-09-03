# v0.7 Release Notes

## Theme

**From one electronic-state tracking path to a branched gauge graph.**

## New modules

- `gaussian_dynamics/gauge_graph.py`
- `gaussian_dynamics/graph_electronic.py`
- `gaussian_dynamics/graph_gaussian.py`
- `gaussian_dynamics/pyscf_gauge_graph.py`

## New theory/docs

- `V07_THEORY.md`
- `V07_DERIVATIONS.md`
- `V07_PYSCF_GRAPH.md`
- `V07_VALIDATION.md`
- `V07_BUILD_VALIDATION.md`

## New examples

- `examples/19_graph_holonomy_ci.py`
- `examples/20_graph_gaussian_gauge_invariance.py`
- `examples/21_pyscf_tbf_centroid_graph.py`

## Main scientific additions

1. Electronic overlap matrices are treated as graph edges rather than a single
   sequential chain.
2. Polar links provide a discrete unitary connection.
3. Wilson loops preserve Berry/non-Abelian holonomy.
4. A spanning-tree gauge isolates loop information on graph chords.
5. Multi-start synchronization smooths noisy local gauges without forcing physical
   holonomy to zero.
6. Electronic operators are transported covariantly, while derivative couplings are
   correctly recognized as connections with an inhomogeneous gauge transformation.
7. TBF electronic vectors are transported to common pair-centroid frames before
   Gaussian pair matrix elements are formed.
