# v0.23.3 algorithmic complexity

Let `R` be replay geometries, `s` retained electronic states, `m` symmetry
manifolds, `d` nuclear coordinates, `F` installed distribution files, and `B`
their total bytes.

## Finite-manifold transport

Each overlap requires an `s x s` SVD, so analysis and polar transport cost
`O(s^3)` time and `O(s^2)` memory. A dense all-pairs replay costs
`O(R^2 s^3)` time and stores `O(R^2 s^2)` raw overlaps and transports. Singular
values add `O(R^2 s)` storage.

## Complete-manifold audit

Projector eigendecompositions and endpoint transport remain `O(s^3)`. Block
retention/leakage comparisons are at most `O(m^2 s^3)` in a direct dense
implementation, with `m <= s`; the electronic-structure calculation remains the
dominant production cost.

## Replay and identity

Canonical manifest hashing is linear in serialized metadata. NPZ validation and
hashing are linear in stored array bytes; validation additionally recomputes the
`O(R^2 s^3)` transport decomposition. Provider/convention fingerprints are linear
in their small canonical JSON records.

## Runtime evidence

Metadata validation is `O(F)`. RECORD-content verification is `O(B)` time with
bounded streaming hash memory and `O(F)` manifest storage. It is a release-evidence
operation, not a per-timestep cost.
