# v0.29.4 — Pre-ML spawning contract

## Purpose

v0.29.4 freezes the boundary between the deterministic Gaussian spawning machinery
and a future machine-learning candidate ranker. It is intentionally a **contract
iteration**, not an ML iteration.

The scientific architecture is:

```text
candidate generator
        |
        v
invariant feature record  ---> future ML ranker (proposal/ranking only)
        |                              |
        |                              v
        +----------------------> selected candidates
                                       |
                                       v
                              exact physics oracle
                                       |
                      residual gain + novelty + metric/
                    conditioning + projection + energy gates
                                       |
                                  spawn/reject
```

The physics oracle remains authoritative. A future model may reduce how many
expensive candidates are evaluated, but it may not bypass physical admission gates.

## Frozen data schema

Each `SpawningRecord` binds four classes of information:

1. **Stable identity** — family, trajectory, event, packet, parent, step and time.
2. **Invariant features** — populations, energy gaps, force norms, NAC Frobenius
   norm, SOC singular values, Berry-curvature norm, Gram conditioning, packet
   overlap, nuclear speed/kinetic energy, residual norm, novelty, projection-loss
   estimate and time step.
3. **Exact oracle supervision** — spawn/no-spawn action, target state, raw and
   normalized residual gain, threshold, reason code, oracle version and final
   admission result.
4. **Candidate provenance** — deterministic candidate kind and serial.

Raw electronic phases, adiabatic eigenvectors, NAC vectors, SOC matrix entries and
other gauge-dependent coordinates are deliberately excluded from the ML boundary.

## Determinism and evidence

All records and manifests use sorted compact JSON with `allow_nan=False` and SHA-256
fingerprints. Dataset hashes are independent of record ordering. Evidence booleans
must be native Python `bool` values, preventing the `numpy.bool_` serialization
ambiguity previously encountered in the curved-TDVP evidence ledger.

A dataset manifest stores:

- schema version;
- split name and record count;
- sorted record fingerprints;
- oracle-version set;
- family/trajectory leakage groups;
- generation commit;
- frozen/unfrozen state;
- deterministic dataset fingerprint.

## Leakage policy

Random candidate-level splitting is forbidden. The minimum grouping unit is
`(family_id, trajectory_id)`. Any group appearing in two different splits causes
`assert_disjoint_manifests` to fail closed.

This matters because many neighboring candidate evaluations share essentially the
same electronic/nuclear state. Candidate-level random splitting would therefore
inflate apparent generalization.

## Invariance bridge

`hermitian_operator_invariants` converts Hermitian electronic operators into trace,
Frobenius norm, eigenvalues and adjacent spectral gaps. These quantities are
unchanged by unitary similarity transformations and provide a small independent
proof that the ML boundary can consume gauge-invariant electronic summaries rather
than basis-dependent matrix entries.

## Evaluation contract

A future ranker must report, at minimum:

- precision;
- recall;
- false-negative rate;
- Brier score / probability calibration;
- mean error in predicted normalized residual gain;
- trajectory-level spawning-event time MAE.

The flat-record evaluator in this iteration implements all metrics except event-time
MAE. Event-time matching is intentionally reserved for the future sequence-level
evaluator, where candidate events can be aligned along complete trajectories.

False-negative rate is explicitly first-class because missing a physically important
spawn can be more damaging than evaluating an extra candidate exactly.

## Non-goals

v0.29.4 does **not**:

- train or select an ML model;
- optimize thresholds on the test split;
- replace the deterministic spawning oracle;
- weaken novelty, metric-conditioning, projection-loss or energy-jump gates;
- claim integration with the unpublished v0.29.2/v0.29.3 curved implementation;
- claim molecular-ab-initio spawning accuracy.

## Integration point

When the curved spawning implementation is moved onto the public development tree,
its candidate-evaluation routine should emit `SpawningRecord` objects immediately
before and after oracle admission. Frozen manifests should then be generated from
large deterministic trajectory campaigns. Only after those datasets pass the
fingerprint, invariance and leakage gates should v0.30/v0.31 introduce a surrogate
ranker.
