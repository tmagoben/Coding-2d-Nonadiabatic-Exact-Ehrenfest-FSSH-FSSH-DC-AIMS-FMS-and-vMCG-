# v0.10 Benchmark Campaign Guide

## A. Exact discretization surface

```python
from gaussian_dynamics import (
    CIPassageConfig,
    run_exact_grid_timestep_surface,
)

rows = run_exact_grid_timestep_surface(
    CIPassageConfig(),
    grid_values=(40, 56, 72),
    dt_values=(0.010, 0.005, 0.0025),
)
```

Inspect changes along **both** axes.

Do not simply choose the finest row and call it converged.

---

## B. Managed Gaussian convergence surface

```python
from gaussian_dynamics import run_managed_parameter_surface

rows = run_managed_parameter_surface(
    dts=(0.010, 0.005, 0.0025),
    spa_orders=(0, 1),
    spawn_action_thresholds=(4e-4, 2e-4, 1e-4),
    max_basis_values=(2, 4, 6),
    overlap_blocks=(0.90, 0.999, 0.9999),
)
```

This is intentionally a multidimensional surface.

Changing several parameters simultaneously does not tell you which approximation
changed the answer.

---

## C. Reduced-density comparison

For the analytic CI model, prefer the global diabatic basis:

```python
from gaussian_dynamics import compare_managed_exact_diabatic_density

comparison = compare_managed_exact_diabatic_density()
```

For a backend without a supplied global diabatic basis, use the common graph frame:

```python
from gaussian_dynamics import compare_managed_exact_common_frame
```

Important outputs:

```text
rho_exact
rho_managed
populations_exact
populations_managed
density_frobenius_error
purity_exact
purity_managed
linear_entropy_exact
linear_entropy_managed
```

This is the preferred v0.10 exact/Gaussian electronic comparison.

---

## D. Wigner initial-condition ensemble

```python
from gaussian_dynamics import run_ci_initial_condition_ensemble

ensemble = run_ci_initial_condition_ensemble(
    nsamples=16,
    seed=12345,
)
```

Increase `nsamples` until the reported SEM is small relative to the physical
difference you are trying to resolve.

---

## E. Interpreting a failure

Suppose:

```text
norm error                   small
exact-grid refinement error  small
Gaussian dt error            small
SPA0-SPA1 difference         small
basis/spawn sensitivity      large
exact-reference error        large
```

The correct conclusion is:

> time propagation and exact reference appear controlled, but the adaptive Gaussian
> representation is not basis converged.

It is **not**:

> the Gaussian method is validated because energy/norm are conserved.

---

## F. Recommended campaign order

Run these in order:

1. exact grid-domain test;
2. exact grid-size refinement;
3. exact timestep refinement;
4. managed timestep refinement with fixed basis controls;
5. SPA0/SPA1 sensitivity;
6. spawning-action sensitivity;
7. overlap-block sensitivity;
8. max-basis sensitivity;
9. pruning/conditioning audit;
10. common-reference density comparison;
11. Wigner ensemble convergence;
12. only then replace the analytic backend with PySCF.

That order prevents electronic-structure uncertainty from masking a nuclear-dynamics
implementation problem.
