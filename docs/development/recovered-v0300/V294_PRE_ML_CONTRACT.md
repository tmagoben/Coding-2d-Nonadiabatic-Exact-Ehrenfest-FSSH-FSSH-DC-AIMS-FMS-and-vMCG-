# v0.29.4 Pre-ML Scientific and Data Contract Hardening

## Scope
v0.29.4 adds **no new physical approximation**. It freezes the interface between the exact v0.29.3 residual-driven curved spawning oracle and any future machine-learning ranker.

The exact evaluator remains authoritative. ML is allowed to propose or rank candidates later; it does not redefine admission physics.

## Frozen oracle target
For candidate `c`, the inherited exact oracle evaluates the full-tangent Schur-complement gain

\[
S_c=C-X^TM^+X,\qquad y_c=q_c-X^TM^+g,\qquad
\Delta_c=y_c^T S_c^+ y_c.
\]

The primary supervised target is the continuous quantity `predicted_residual_reduction = Delta_c`, with `relative_residual_reduction`, Schur rank/condition and exact rejection reasons retained as auxiliary targets.

## ML boundary

```text
trajectory state
   -> deterministic candidate generator
   -> exact v0.29.3 candidate evaluator
   -> SpawnOracleRecordV294
   -> versioned JSONL dataset
   -> future ML ranker
   -> exact physics admission (retained)
```

Candidate generation and candidate evaluation are deliberately separate. A future ML model may replace or accelerate ranking without changing the exact evaluator.

## Canonical state contract
Every packet is serialized as

\[
(q_I,p_I,\Gamma_I,B_I,c_I,\mathrm{id})
\]

with packet rows sorted by stable ID. The state fingerprint is provenance, not an ML feature. Path-family identity is also recorded in state provenance.

## Invariant feature contract
The first feature schema deliberately excludes raw electronic coefficients and raw gauge-dependent residual vectors. It records scalar/canonical quantities including:

- source packet population;
- width-metric displacement \(\Delta q^T\Gamma\Delta q\);
- inverse-width momentum displacement \(\Delta p^T\Gamma^{-1}\Delta p\);
- Euclidean displacement norms;
- sorted width and chirp eigenvalues;
- width log determinant;
- chirp Frobenius norm;
- TDVP metric rank and condition;
- TDSE residual norm;
- displacement kind.

These are intended to be invariant to packet ordering and electronic gauge choices and substantially safer under orthogonal coordinate transformations than raw component features.

## Rejection taxonomy
Exact reasons are retained rather than collapsing all failures into one binary label. The vocabulary includes low residual gain, no new tangent direction, ill-conditioning, projection loss, energy jump, duplicate/non-novel candidates, degenerate width eigenspaces, path/transport failures and numerical failures.

## Leakage policy
Train/validation/test assignment is done by **trajectory group ID**, never by individual candidate row. Every candidate emitted by one trajectory therefore lands in the same split. Later molecular datasets should use still harder group keys (molecule/model family/trajectory seed) where appropriate.

## Determinism
Runtime cost is measured separately from immutable scientific records. Wall-clock timing is not hashed into dataset labels. Two recomputations of the fixture produce byte-identical JSONL and the same dataset/evidence fingerprint.

## Future ML safety boundary
Initial ML integration must remain

\[
\boxed{\text{ML shortlist/rank} \rightarrow \text{exact v0.29.3/v0.30 admission}}
\]

A deterministic sample of ML-rejected candidates is reserved for exact false-negative audits and distribution-shift detection.
