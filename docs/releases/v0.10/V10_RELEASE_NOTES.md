# v0.10 Release Notes

Version 0.10 is the **benchmark-campaign and observable-correctness release**.

New modules:

```text
initial_conditions.py
benchmark_metrics.py
error_budget.py
benchmark_acceptance.py
benchmark_campaign.py
ensemble_benchmark.py
campaign_io.py
electronic_observables.py
reference_comparison.py
```

Major scientific changes:

1. exact and Gaussian calculations can be compared through the same fixed-frame
   reduced electronic density matrix;
2. reduced-state purity and entropy expose missing electron-nuclear branching;
3. Gaussian Wigner initial conditions are sampled reproducibly;
4. exact `grid x dt` surfaces are first-class objects;
5. managed `dt x SPA x spawning x basis x overlap-block` surfaces are first-class
   objects;
6. sensitivity/error budgets are reported without pretending correlated errors are
   additive;
7. benchmark acceptance criteria are explicit;
8. repeated spawning is available as an opt-in extension while preserving earlier
   API behavior.

The central philosophy is simple:

> v0.10 is allowed to demonstrate that the present Gaussian approximation is not yet
> converged for a demanding CI passage.

That information determines what v0.11 should improve.


## Build result

The cumulative automated suite for this release is:

```text
111 passed
```

The compact release benchmark intentionally fails the configured exact-reference
population criterion.

For the default near-CI passage:

```text
exact global-diabatic reduced populations:
[0.22600611, 0.77399389]

managed reference reduced populations:
[0.95027840, 0.04972160]

exact reduced-state purity:
0.67620820

managed reduced-state purity:
0.97488525
```

The compact sensitivity budget identifies **basis size** as the largest controlled
proxy among the tested refinement axes, while the total exact discrepancy is still
larger.  This is evidence that the current small-basis spawning representation is not
yet adequate for the stronger passage benchmark.
