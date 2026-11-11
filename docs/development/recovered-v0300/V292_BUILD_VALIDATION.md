# v0.29.2 Build Validation

## Classification
Development checkpoint layered on the certified v0.27.0 substrate and the retained v0.29.1 mathematical contract. It is not a public release and package metadata remains unchanged.

## Implemented surface
- constant anti-Hermitian non-Abelian connection and straight Wilson candidate sections;
- complete curved-state finite-difference TDVP Jacobian over (Re c, Im c, q, p, E, B);
- full-tangent Schur-complement residual-gain score;
- full-SVD curved wavefunction projection;
- curved merge, prune, spawn, and exact no-op lifecycle ordering;
- curved packet-subspace overlap for merge admission;
- projection-loss and energy-jump gates;
- coefficient-active / shape-dormant newborn state;
- population + condition + velocity-amplification activation machinery;
- deterministic candidate ordering and degenerate-width fail-closed inheritance.

## Focused validation
- New v0.29.2 tests: **9/9 passed**.
- Replayed inherited v0.27 correlated TDVP/basis/validation tests: **22/22 passed**.
- Focused total: **31/31 passed**.
- Python compileall over package/examples/tests: passed.

## Independent residual-gain check
Predicted Schur-complement reduction: `1.1783841637889198e-02`

Independent augmented real least-squares reduction: `1.1783841637889186e-02`

Relative disagreement: `1.030e-15`

## Spawn projection
- relative projection loss: `4.361e-31`
- fidelity: `1.0000000000000000`
- newborn coefficient norm: `4.024e-16`
- exact no-event object identity: `True`

## Evidence
- evidence fingerprint: `40f9fa794a2f1e8991871291e485d798bddc7a7b1178d79c1ea2e73ce54bae99`
- evidence JSON SHA-256: `4032d8d4824ca5f4a43a956ac5cd1e6fb1a3b36f2ecd922c42e0a8b9a8cba5fc`

## Honest boundaries
This checkpoint does not claim coordinate-dependent connections, arbitrary path families, bundle patch transitions, live molecular derivative couplings/SOC, multi-event AIMS branching, ML spawning, or nuclear-spin dynamics. A complete historical v0.29.0/v0.29.1 source replay is not claimed because those development source files are retained outside the certified public v0.27 tree; their retained evidence records remain the scientific baseline.
