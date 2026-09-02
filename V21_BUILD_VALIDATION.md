# v0.21 Build Validation

Validated on 2026-08-13.

## Source validation

```text
351 Python files parsed successfully with Python AST.
```

## Cumulative regression suite

```text
240 passed in 35.09 s
```

The suite remains cumulative across the earlier Gaussian-dynamics releases.

New v0.21 coverage includes:

```text
general complex electronic-operator contract
physical Hamiltonian-derivative matrices
complex anti-Hermitian derivative connections
smooth coordinate-dependent U(s) gauges
gauge-invariant generalized force
cross-geometry complex-overlap covariance
full-subspace Procrustes alignment
Wilson-loop spectrum invariance
block S/H/T complex-gauge covariance
gauge-invariant block sparse edge scoring
all-edge sparse block equality with dense block reference
time-dependent gauge-equivalent propagation
second-order gauge-equivalence convergence
arbitrary five-state block construction
2/4/8-state sparse block scaling
dynamic block graph edge entry and exit
complete v0.21 release acceptance
```

## Canonical release campaign

Machine-readable output:

```text
results/v021_complex_block_framework_campaign.json
```

Acceptance:

```text
passed = True
```

Every configured release check passes.

## Complex point covariance

```text
H relative error:
0.0

maximum physical-Hamiltonian-derivative error:
0.0

gauge-invariant force error:
8.673617379884035e-19
```

## Complex block covariance

```text
S relative error:
1.7399964775418827e-16

H relative error:
1.7933278321452532e-16

T relative error:
1.7899468791987713e-16

maximum block-edge score change:
5.551115123125783e-17
```

The \(T\) test includes the inhomogeneous term

$$
\mathcal G^\dagger S\dot{\mathcal G},
$$

so this is a genuinely coordinate-dependent gauge test.

## Gauge-equivalent propagation

```text
dt=0.020:
error = 2.047450337289394e-08

dt=0.010:
error = 5.118673197623611e-09

dt=0.005:
error = 1.2796717006233196e-09
```

Observed orders:

```text
[1.999986653128339, 1.9999961654805536]
```

Minimum observed order:

```text
1.999986653128339
```

The measured behavior is essentially second order.

## Degenerate-subspace validation

```text
tracked subspace dimension:
8

minimum singular value:
0.9999999999999998

anti-Hermitian residual after Procrustes alignment:
1.967477117880543e-15
```

## Wilson-loop validation

```text
maximum complex-gauge eigenphase change:
9.43689570931383e-16
```

## Block sparse convergence

Zero local omission budget:

```text
S error:
1.3757894522835919e-07

H error:
1.5178432514327997e-06

T error:
1.0997753094998538e-06
```

Threshold convergence:

```text
{'H_error': True, 'S_error': True, 'T_error': True}
```

Budget convergence:

```text
{'H_error': True, 'S_error': True, 'T_error': True}
```

## Dynamic topology

```text
entered edges:
15

exited edges:
9

maximum active edges:
15

total exact pair checks:
749
```

## Electronic-state dimension

```text
tested state dimensions:
[2, 4, 8]

H_nnz / s^2:
[114.0, 114.0, 114.0]

relative scaling error:
0.0
```

## PySCF and spin scope

```text
PySCF installed in build environment:
False

real PySCF v0.21 trajectory validated:
False

SOC Hamiltonian introduced:
False
```

v0.21 validates the **general complex framework**, not physical spin-orbit dynamics.

## Scientific limitations

The release does not claim:

```text
production AIMS
complete production molecular matrix elements
real PySCF v0.21 runtime validation
physical spin-orbit coupling dynamics
production asynchronous electronic-structure scheduling
```
