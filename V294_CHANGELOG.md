# v0.29.4 changelog

## Added

- Immutable `SpawningRecord` schema for future ML-assisted spawning datasets.
- Gauge-invariant electronic/nuclear feature contract.
- Exact deterministic spawning-oracle supervision contract.
- Canonical JSON serialization and SHA-256 record/dataset fingerprints.
- Frozen dataset manifests with family/trajectory leakage groups.
- Fail-closed split-leakage auditing.
- Unitary-similarity invariant Hermitian-operator feature bridge.
- Prediction evaluator emphasizing recall, false-negative rate, calibration and
  normalized-residual-gain error.
- Focused tests, build-validation record, acceptance gates and pre-ML roadmap.

## Scientific boundary

v0.29.4 is a standalone pre-ML contract layer. It does not train a model, replace the
physics spawning oracle, or claim integration with unpublished v0.29.2/v0.29.3 curved
dynamics. The intended architecture is strictly:

```text
ML proposes/ranks -> exact physics evaluates -> physical gates admit/reject
```
