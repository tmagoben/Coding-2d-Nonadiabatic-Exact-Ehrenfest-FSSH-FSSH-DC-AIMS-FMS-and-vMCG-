# v0.29.4 Build Validation

## Classification
Development checkpoint. No new physics claim and no replacement for the certified v0.27 release.

## Added
- `gaussian_dynamics/pre_ml_contract_v294.py`
- `gaussian_dynamics/v294_benchmark.py`
- `tests/test_v0294_pre_ml_contract.py`
- `results/v0294_pre_ml_evidence.json`
- `results/v0294_spawn_dataset_v1.jsonl`
- `V294_PRE_ML_CONTRACT.md`

## Acceptance campaign
Focused tests are executed in three independently reported groups:

- v0.29.4 pre-ML contract: **7/7 passed**;
- v0.29.2-v0.29.3 curved spawning/geometry: **19/19 passed**;
- inherited v0.27 correlated TDVP/basis/validation: **22/22 passed**.

Total focused checkpoint: **48/48 passed**.

## Deterministic evidence
`python -m gaussian_dynamics.v294_benchmark` writes a three-record schema fixture spanning train, validation and test trajectory groups.

The deterministic dataset is byte-identical across repeated recomputation. Timing is recorded only as non-hashed performance metadata.

Current evidence and dataset fingerprints are stored in `results/v0294_pre_ml_evidence.json`.

## Explicit boundaries
- This fixture is a schema/contract demonstration, not a production training corpus.
- No ML model is trained in v0.29.4.
- No ML prediction is permitted to override exact spawning admission.
- No live molecular NAC/SOC trajectory is claimed here.
- Raw gauge-dependent electronic coefficient components are deliberately excluded from the initial ML feature schema.
