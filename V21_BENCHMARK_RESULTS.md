# v0.21 Benchmark Results

Machine-readable output:

```text
results/v021_complex_block_framework_campaign.json
```

## Complex point covariance

```text
H relative error:
0.0

maximum dH relative error:
0.0

force error:
8.673617379884035e-19
```

## Full block covariance

```text
S:
1.7399964775418827e-16

H:
1.7933278321452532e-16

T:
1.7899468791987713e-16

maximum edge-score error:
5.551115123125783e-17
```

## Time-dependent gauge propagation

```text
dt=0.02  error=2.047450337289394e-08  norm_drift=6.570369803782228e-09
```
```text
dt=0.01  error=5.118673197623611e-09  norm_drift=1.6426023874416273e-09
```
```text
dt=0.005  error=1.2796717006233196e-09  norm_drift=4.10651956883612e-10
```

Observed orders:

```text
[1.999986653128339, 1.9999961654805536]
```

## Full-subspace tracking

```text
states:
8

minimum singular value:
0.9999999999999998

anti-Hermitian alignment residual:
1.967477117880543e-15
```

## Wilson loop

```text
maximum eigenphase difference:
9.43689570931383e-16
```

## Dynamic sparse topology

```text
entered edges:
15

exited edges:
9

maximum active:
15

final active:
6
```

## Electronic-state scaling

```text
s=2  dimension=48  H_nnz=456  density=0.19791666666666666
```
```text
s=4  dimension=96  H_nnz=1824  density=0.19791666666666666
```
```text
s=8  dimension=192  H_nnz=7296  density=0.19791666666666666
```

`H_nnz / s^2`:

```text
[114.0, 114.0, 114.0]
```

## Acceptance

```text
passed = True
```

All configured v0.21 checks pass.
