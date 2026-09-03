# v0.21 Release Notes

v0.21 is the **complex/block/gauge-ready generalization release**.

It deliberately introduces **no spin physics**.

## New modules

```text
electronic_operator_v21.py
complex_gauge_v21.py
subspace_tracking_v21.py
wilson_loop_v21.py
synthetic_operator_provider_v21.py
block_sparse_molecular_v21.py
block_dynamics_v21.py
v21_benchmark.py
```

## Major changes

```text
general complex electronic operator contract
arbitrary electronic block dimension
complex U(s) gauge covariance
physical Hamiltonian-derivative operators
gauge-invariant block sparse edge score
full-subspace Procrustes alignment
Wilson-loop spectrum validation
block-sparse molecular S/H/T
time-dependent complex-gauge propagation validation
dynamic sparse topology stress
2/4/8-state regression coverage
curated derivation/complexity/validation documentation
```

## Key result

The same prescribed molecular Gaussian dynamics in two coordinate-dependent complex
gauges converges together at essentially second order:

```text
observed orders:
[1.999986653128339, 1.9999961654805536]

finest gauge-mapped coefficient error:
1.2796717006233196e-09
```

## Scope

The framework can accept a future SOC Hamiltonian through the electronic operator
backend, but v0.21 itself remains fully usable for ordinary spin-free nonadiabatic
dynamics.

No production AIMS or real PySCF v0.21 trajectory is claimed.


## Cumulative validation

```text
351 Python files parsed
240 tests passed
```

See `V21_BUILD_VALIDATION.md` for the complete release contract.
