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

## Documentation and mathematics

For Markdown rendered on GitHub:

- use `$...# Contributing

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

 for inline mathematics and `$...$` for display mathematics;
- prefer `\\mathrm{...}` to the legacy `\\rm` declaration;
- keep mathematical notation outside code spans unless the literal source is being discussed;
- preserve the scientific meaning and historical claim boundary of versioned release documents.

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
