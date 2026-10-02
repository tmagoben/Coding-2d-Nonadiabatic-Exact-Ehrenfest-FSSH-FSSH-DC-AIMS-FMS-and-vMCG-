# v0.29.5 — Curved-oracle integration and dataset certification

## Scope

This iteration turns the v0.29.4 pre-ML schema into a campaign-level certification
layer. It still does **not** train a surrogate and it does not allow learned output to
admit a Gaussian spawn.

## Implemented

- Deterministic train/validation/test/audit assignment at the complete
  `family_id + trajectory_id` group level.
- Exact replay receipts that compare complete invariant features, candidate
  provenance and oracle decision payloads.
- Invariance receipts for electronic-gauge, packet-permutation and coordinate
  transformations supplied by the physics implementation.
- Campaign manifests bound to the generation commit and v0.29.4 record hashes.
- A campaign fingerprint that changes if manifests, replay evidence, invariance
  evidence or pass/fail state changes.
- Canonical JSONL serialization for downstream dataset materialization.

## Fail-closed rules

A campaign is not certified when:

1. the physics oracle cannot exactly replay its stored decision;
2. replay mutates the invariant feature record;
3. a registered physical symmetry changes the ML-facing record;
4. family/trajectory leakage occurs across splits;
5. provenance is missing or non-deterministic.

## Integration boundary

The public repository still does not contain the advanced curved v0.29.2/v0.29.3
TDVP/spawning engine. Therefore this branch implements the **integration harness**,
not a fabricated curved oracle.

When that engine is recovered into the repository, three adapters are required:

1. emit a v0.29.4 `SpawningRecord` for every exact candidate evaluation;
2. expose an exact replay callback that recomputes the candidate and oracle decision;
3. expose physical gauge/permutation/coordinate transforms whose re-evaluated
   records are passed through `audit_invariant_record_transform`.

Only a campaign for which `certify_campaign(...).passed is True` is eligible to
be frozen for future surrogate training.

## Next acceptance step

Run the combined v0.29.4 + v0.29.5 focused test suite, then connect the recovered
curved spawning oracle and require real (not identity-placeholder) invariance
transforms before calling v0.29.5 scientifically integrated.
