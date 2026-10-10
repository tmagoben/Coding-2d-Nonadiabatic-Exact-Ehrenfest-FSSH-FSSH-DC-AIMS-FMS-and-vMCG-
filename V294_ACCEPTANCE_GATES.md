# v0.29.4 acceptance gates

v0.29.4 is accepted only when all gates below are satisfied.

## Contract gates

- [x] Schema version is frozen at `0.29.4`.
- [x] Candidate records contain stable family/trajectory/event/packet identity.
- [x] ML-facing electronic information is reduced to gauge-invariant summaries.
- [x] Exact oracle labels retain residual gain, threshold, reason and admission state.
- [x] Record and dataset serialization rejects NaN and non-finite values.
- [x] Record fingerprints are deterministic.
- [x] Dataset fingerprints are independent of record ordering.
- [x] Train/validation/test/audit splitting fails closed on family/trajectory leakage.
- [x] Future prediction evaluation makes false-negative rate and calibration explicit.
- [x] An ML surrogate has no spawn-admission authority.

## Validation gates

- [x] 14/14 focused tests pass locally.
- [x] Python bytecode compilation succeeds for module and tests.
- [x] Unitary electronic-basis changes preserve the tested Hermitian invariants.
- [x] `numpy.bool_` is rejected at the evidence boundary.
- [x] Non-Hermitian operator input is rejected.

## Integration gates intentionally open

These are **not** prerequisites for merging the standalone contract, but they must be
closed before a learned ranker is used in production dynamics:

- [ ] Publish/integrate the v0.29.2/v0.29.3 curved spawning oracle on the development tree.
- [ ] Emit `SpawningRecord` at every candidate evaluation in deterministic trajectory campaigns.
- [ ] Freeze leakage-safe train/validation/test manifests from multiple molecular/model families.
- [ ] Add sequence-level spawn-event matching and `event_time_mae` evaluation.
- [ ] Demonstrate gauge/permutation/path invariance of generated datasets end to end.
- [ ] Benchmark candidate-evaluation savings with the exact physics oracle still downstream.
- [ ] Establish uncertainty/abstention rules before any learned ranking can suppress exact evaluations.

The architecture remains `ML proposes/ranks -> physics verifies`.
