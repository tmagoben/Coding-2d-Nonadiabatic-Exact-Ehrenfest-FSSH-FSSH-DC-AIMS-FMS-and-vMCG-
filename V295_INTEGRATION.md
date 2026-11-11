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

At the original v0.29.5 integration, the public repository did not contain the
curved v0.29.2/v0.29.3 TDVP/spawning engine. The v0.30.1 development branch restores
that research implementation from the verified v0.30.0 checkpoint. This file still
describes the **integration harness**: restoring the engine does not certify it.

Three adapters remain required. In particular, the recovered package-level
`SpawnOracleRecordV294` and the top-level `SpawningRecord` are different contracts:

1. emit a v0.29.4 `SpawningRecord` for every exact candidate evaluation;
2. expose an exact replay callback that recomputes the candidate and oracle decision;
3. expose physical gauge/permutation/coordinate transforms whose re-evaluated
   records are passed through `audit_invariant_record_transform`.

Only a campaign for which `certify_campaign(...).passed is True` is eligible to
be frozen for future surrogate training.

## Next acceptance step

Connect the recovered curved spawning oracle through an explicit schema adapter
and require physical invariance transforms before calling v0.29.5 scientifically
integrated. Passing the focused test suites verifies their separate interfaces;
it is not sufficient evidence for an end-to-end curved dataset campaign.
