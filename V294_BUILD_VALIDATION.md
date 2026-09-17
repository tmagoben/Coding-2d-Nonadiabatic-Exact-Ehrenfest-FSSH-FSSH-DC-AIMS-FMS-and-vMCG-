# v0.29.4 build validation

## Scope

This evidence record covers the standalone pre-ML spawning contract on branch
`feature/v0294-preml-contract`. It does not certify integration with unpublished
v0.29.2/v0.29.3 curved dynamics code.

## Local validation

The exact module and test snapshot was exercised with:

```bash
python -m pytest -q
python -m py_compile v294_preml_contract.py test_v294_preml_contract.py
```

Result:

```text
14 passed in 0.07s
```

The tests cover:

1. frozen schema identity (`0.29.4`);
2. canonical JSON independent of mapping insertion order;
3. deterministic and physics-sensitive record fingerprints;
4. rejection of NaN/non-finite features;
5. normalized electronic populations;
6. spawn/no-spawn target-state consistency;
7. native-Python boolean evidence boundaries;
8. order-independent dataset fingerprints;
9. fail-closed family/trajectory leakage detection;
10. valid disjoint train/test groups;
11. unitary-similarity invariance of Hermitian electronic summaries;
12. rejection of non-Hermitian electronic operators;
13. false-negative and residual-score evaluation metrics;
14. JSON-serializable, self-consistent manifests.

## Local file fingerprints

Before upload, the tested files had SHA-256 values:

```text
0eecd432f3dc8176b060bf6c34b651327a1ca6a6735c385a4b3cc23519c46931  v294_preml_contract.py
6ba7c34833faa76354a629dee3cda0eb39a08ac2d956e2f3358f1e0a0621042b  test_v294_preml_contract.py
```

The design document was validated as documentation only and is not part of the
runtime test fingerprint.

## Scientific claim boundary

Passing these tests establishes a deterministic data/provenance/evaluation boundary
for future spawning surrogates. It does **not** establish surrogate accuracy,
curved-TDVP integration, molecular accuracy, or authority for an ML model to admit a
spawn. Exact physical admission remains downstream of any future learned ranking.
