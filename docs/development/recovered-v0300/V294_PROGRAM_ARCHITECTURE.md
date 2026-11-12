# v0.29.4 Program Architecture

```text
v0.29.3 state + connection + path family
        |
        v
v0.27 intrinsic candidate generator
        |
        v
exact full-tangent v0.29.3 residual oracle
        |
        +--> exact decision/rejection taxonomy
        +--> invariant feature extractor
        +--> provenance/config fingerprints
        |
        v
SpawnOracleRecordV294
        |
        v
trajectory-group split + deterministic JSONL
        |
        v
future ML ranker (not implemented)
        |
        v
exact physics verification remains mandatory
```

The data contract is intentionally downstream of the exact scientific oracle. It must not be imported into the TDVP mathematics to make the solver depend on ML conventions.
