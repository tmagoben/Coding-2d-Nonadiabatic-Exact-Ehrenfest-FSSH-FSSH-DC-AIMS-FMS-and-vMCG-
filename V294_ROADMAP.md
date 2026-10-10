# v0.29.4 to v0.31 roadmap

## v0.29.4 — contract and evidence boundary

Freeze deterministic records, invariant features, exact-oracle supervision,
leakage-safe manifests, fingerprints and evaluation requirements. No model training.

## v0.30.x — dataset generation on curved molecular dynamics

After the curved spawning implementation is public and integrated:

- emit one immutable record per candidate evaluation;
- cover multiple path geometries, gauge choices and molecular/model families;
- record exact residual improvement and every downstream admission gate;
- freeze trajectory/family-disjoint train, validation and test manifests;
- add sequence-level event matching and dataset invariance audits.

The dataset target should remain continuous normalized residual improvement whenever
possible, with spawn/no-spawn retained as a derived/admission label.

## v0.31 — ML-assisted candidate ranking

Introduce the first deliberately small, interpretable surrogate. Compare a tree-based
regressor/ranker and a compact MLP only after the dataset is frozen. The model may
rank candidates or abstain; it may not override the exact physical oracle.

Primary evaluation questions:

1. How many exact candidate evaluations can be avoided at fixed oracle recall?
2. What is the false-negative rate for physically admitted spawns?
3. Is predicted residual gain calibrated across molecular families and curvature regimes?
4. Does ranking preserve gauge, packet-permutation and path invariance?
5. Does an out-of-distribution/uncertain case reliably fall back to exact evaluation?

A direct learned `spawn/no_spawn` authority remains out of scope until these questions
are answered with independent evidence.
