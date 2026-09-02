# v0.18 Release Notes

v0.18 is the **dynamics and convergence-completeness release**.

## New modules

```text
wavefunction_metrics_v18.py
defect_candidates_v18.py
sampled_sparse_audit_v18.py
convergence_complexity_v18.py
convergence_complete_dynamics_v18.py
convergence_campaign_v18.py
convergence_worker_v18.py
v18_benchmark.py
```

## Major additions

```text
full spinor Gaussian wavefunction reconstruction
global-phase-aligned fidelity and L2 error
nuclear-density L2 and total variation
centroid and covariance errors
exact projected-reference separation
physical-time normalized adaptive controls
batched candidate-grid residual ranking
sampled omitted-edge sparse audits
initial/final dense sentinels only
basis-size convergence ladder
timestep self-convergence
sparse-edge-budget convergence
adaptive-growth threshold sensitivity
fresh-process convergence workers
```

## Canonical numerical result

```text
projected-reference fidelity:
0.982566093411826

phase-aligned L2:
0.13232747836123407

nuclear-density L2:
0.052341235444456596

projected reduced-density error:
0.00010573932284646514

norm drift:
1.2434515149761793e-06
```

## Second-order timestep convergence

```text
observed order:
1.9974143869640382
```

## Basis improvement

```text
Nmax 10 error:
0.19119002020634157

Nmax 13 error:
0.13232747836123407

relative improvement:
30.79 %
```

## Audit and memory improvements

```text
dense audit pair-work reduction vs v0.17:
71.15 %

candidate-grid peak-memory reduction:
97.55 %
```

## Release acceptance

```text
passed = True
```

The cumulative automated regression suite reports:

```text
212 passed
```

Full details are recorded in `V18_BUILD_VALIDATION.md`.

## Scope

v0.18 remains a controlled analytic-LVC research prototype.

It does not yet claim production AIMS, production molecular PySCF dynamics, or
spin-orbit-coupled dynamics.
