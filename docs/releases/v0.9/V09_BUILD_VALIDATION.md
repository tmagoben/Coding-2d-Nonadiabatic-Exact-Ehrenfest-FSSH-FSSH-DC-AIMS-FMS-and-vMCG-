# v0.9 Build Validation Report

Validated on 2026-08-12.

## Compilation

All Python files under `gaussian_dynamics/`, `examples/`, and `tests/` compile successfully with `python -m py_compile`.

## Automated regression suite

```text
........................................................................ [ 73%]
..........................                                               [100%]
98 passed in 2.37s
```

The cumulative suite contains all v0.1-v0.8 tests plus v0.9 tests for:

- exact SPA1 integration of a linear field between multidimensional Gaussians;
- Hermiticity of the graph SPA1 matrix;
- duplicate-basis projection;
- condition-number-controlled TBF pruning;
- canonical orthogonalization;
- timestep-aware integrated coupling exposure;
- synthetic second-order convergence-order recovery;
- exact 2D adiabatic-population normalization;
- managed graph-AIMS spawning and generalized norm conservation;
- pruning an initially singular/redundant TBF basis before propagation.

## Representative v0.9 outputs

### SPA0 versus first-order electronic Taylor layer

```text
SPA0 final populations: [5.21451302e-09 1.00000002e+00]
SPA1 final populations: [1.44300534e-09 1.00000001e+00]
L2 difference: 7.183820373960817e-09

The difference is a controlled approximation diagnostic for this short analytic CI run.
```

### Basis pruning

```text
condition before: 19999998344.19272
condition after:  1.0
removed indices:  (1,)
projection loss:  5.000000413701855e-11
projected C:      [1.+0.j]
```

### Integrated coupling-action spawning

```text
dt=    0.01  trigger_time=    0.05  accumulated_action=     0.1
dt=   0.005  trigger_time=   0.055  accumulated_action=    0.11
dt=  0.0025  trigger_time=    0.05  accumulated_action=     0.1

The integrated |v.d| dt criterion is substantially less tied to an arbitrary per-step threshold.
```

### Managed graph-AIMS versus exact 2D CI reference

```text
Exact final adiabatic populations:   [9.23112112e-07 9.99999077e-01]
Managed graph-AIMS populations:     [5.21451302e-09 1.00000002e+00]
Population L2 error:                 1.312934948476658e-06
Exact norm:                          0.9999999999999973
Managed generalized norm:            1.0000000000000002
Spawn events:                        [{'kind': 'spawn', 'step': 2, 'time': 0.0004, 'parent_uid': 0, 'child_uid': 1, 'target_state': 0, 'integrated_coupling_action': 3.3661986522897427e-06, 'instantaneous_coupling_rate': 0.008415381649645117}]
```

### Timestep refinement

```text
dt values:
[0.0008 0.0004 0.0002]

final populations:
[[5.21439249e-09 1.00000002e+00]
 [5.21451300e-09 1.00000002e+00]
 [5.21451302e-09 1.00000002e+00]]

successive population-vector errors:
[3.81044895e-13 2.88776072e-17]

observed refinement orders:
[13.68772216]
```

## Interpretation of the short exact-reference comparison

For the supplied short-time analytic CI benchmark, the final population-vector difference between the managed Gaussian run and the exact 2D split-operator reference is approximately `1.31e-6`.

This is an illustrative regression case, not a universal AIMS accuracy statement. Longer propagation, stronger branching, and broader basis/refinement studies are required before making physical accuracy claims.

## SPA terminology

The v0.9 `order=0` layer is a centroid/zeroth-order saddle-point electronic approximation.

The `order=1` layer adds the analytically integrated first Taylor term. It is deliberately documented as an **SPA1 electronic Taylor layer**, not as a complete production implementation of every first-order AIMS Hamiltonian term near a conical intersection.

## PySCF status

The explicit PySCF backend, many-electron state tracking, gauge graph, and incremental snapshot graph from v0.5-v0.8 remain included.

PySCF is not installed in this build environment, so those binary-backed examples are not executed here. Their API contracts remain covered by the inherited fake-backend and overlap-engine tests.

See `V09_PYSCF_CONVERGENCE.md` for the recommended electronic-structure convergence sequence before using the managed Gaussian dynamics on a molecule.
