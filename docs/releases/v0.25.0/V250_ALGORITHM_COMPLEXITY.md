# v0.25.0 algorithmic complexity

Let `d` be the number of nuclear coordinates and `s` the retained electronic-state
count. Provider cost is denoted `C_provider`; it remains backend dependent.

## One step

- Two operator snapshots: `2 C_provider`.
- One cross-geometry overlap: backend-dependent `C_overlap`.
- Two dense force contractions over all coordinates: `O(d s^2)` time.
- One constant dense mass solve: `O(d^3)` without cached factorization, or `O(d^2)`
  per step after factoring a fixed mass matrix; `O(d^2)` storage.
- One `s x s` SVD and associated polar diagnostics: `O(s^3)` time and `O(s^2)`
  storage.
- Two Hermitian eigendecompositions/exponentials: `O(s^3)` time and `O(s^2)`
  storage.
- Step-receipt storage: `O(d s^2+s^2+d^2)`.

For `N` retained receipts, trajectory storage is linear in `N`. A streaming runner
could reduce this to constant memory, but v0.25.0 intentionally keeps all receipts
for deterministic validation and tamper detection.

## Numerical choice consequences

Using the SVD is not an avoidable extra if trajectory readiness must be measured:
the same factorization yields both `W=UV^dagger` and every singular-value quality
gate. A raw-overlap multiplication would be cheaper but nonunitary and physically
incorrect for amplitude transport.

The explicit constant-mass check is cheap. The future full TDVP will instead require
assembling and solving a variational metric system, likely with dense worst-case cost
`O(n_z^3)` per nonlinear iteration for `n_z` variational degrees of freedom. That is
why it is a separate implementation milestone rather than an option on this Verlet
runner.
