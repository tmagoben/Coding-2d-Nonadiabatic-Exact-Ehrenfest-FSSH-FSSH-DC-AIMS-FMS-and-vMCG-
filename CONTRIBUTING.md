# Contributing

## Scientific requirements

A numerical change should include at least one of:

- an analytic identity test;
- an exact-grid or dense-reference comparison;
- a gauge/invariance test;
- a convergence test;
- a deterministic regression test.

Do not accept a scientific algorithm change based only on visual agreement.

## Code style

Prefer:

- small functions;
- explicit array shapes;
- explicit physical conventions;
- linear solves instead of matrix inverses;
- deterministic random seeds in tests;
- comments explaining mathematics rather than Python syntax.

## New physics

A new physical model should enter through a narrow interface.

A future SOC implementation should provide a valid complex Hermitian electronic
operator instead of modifying the Gaussian propagator to assume spin everywhere.

## Validation

Before committing:

```bash
python -m pytest -q
```

For release changes, recompute the versioned benchmark campaign and document the
machine-readable result.
