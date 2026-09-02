# v0.18 Validation Contract

v0.18 is accepted only if reduced-state accuracy, full-wavefunction accuracy,
discretization convergence, basis convergence, sparse convergence, and adaptive-growth
sensitivity all pass simultaneously.

## Canonical full-wavefunction requirements

```text
projected-reference fidelity              >= 0.980
phase-aligned wavefunction L2 error       <= 0.140
nuclear-density L2 error                  <= 0.060
centroid L2 error                         <= 0.003
covariance Frobenius error                <= 0.015
```

Measured:

```text
fidelity:
0.982566093411826

wavefunction L2:
0.13232747836123407

density L2:
0.052341235444456596

centroid error:
0.001269391116081437

covariance error:
0.01269760151753039
```

All pass.

## Reduced-state and conservation requirements

```text
projected reduced-density error     <= 2e-4
target reduced-density error        <= 0.035
generalized norm drift              <= 1e-4
condition number                    <= 1e4
```

Measured:

```text
projected reduced-density error:
0.00010573932284646514

target reduced-density error:
0.03329249794783041

norm drift:
1.2434515149761793e-06

condition number:
6509.218903498147
```

The condition ceiling is intentionally looser than v0.17 because the 13-Gaussian basis
improves full-wavefunction fidelity by more than 30% relative to the 10-Gaussian
budget, while remaining numerically stable.

## Exact reference consistency

Initial exact projected-target fidelity:

```text
0.8822514544600696
```

Final exact projected-target fidelity:

```text
0.8822514544600707
```

Maximum drift:

```text
1.6653345369377348e-15
```

Required:

```text
<= 1e-9
```

This validates the representation-error decomposition.

## Trajectory fidelity

Stored projected-reference fidelities:

```text
t=0.0    F=1.000000000000    L2=0.000000000
```
```text
t=0.1    F=0.999092291223    L2=0.030131625
```
```text
t=0.2    F=0.998002493180    L2=0.044704645
```
```text
t=0.3    F=0.996445507167    L2=0.059646098
```
```text
t=0.4    F=0.993684031331    L2=0.079535987
```
```text
t=0.5    F=0.989218592711    L2=0.103974154
```
```text
t=0.6    F=0.982566093412    L2=0.132327478
```

Minimum fidelity:

```text
0.9825660934118264
```

Required:

```text
>= 0.975
```

## Timestep convergence

Successive phase-aligned wavefunction differences:

```text
dt 0.010 -> 0.005:
0.0006101685310847837

dt 0.005 -> 0.0025:
0.00015281576549629837
```

Observed order:

```text
1.9974143869640382
```

Required:

```text
>= 1.5
```

The measured value is essentially second order.

## Basis convergence

Errors:

```text
Nmax=10  L2=0.19119002020634157
```
```text
Nmax=11  L2=0.17750196410133656
```
```text
Nmax=12  L2=0.14492430172716658
```
```text
Nmax=13  L2=0.13232747836123407
```

Requirements:

```text
strictly decreasing error
relative improvement >= 25%
```

Measured relative improvement:

```text
30.79 %
```

## Sparse-edge-budget convergence

```text
B_local=0.03  L2=0.14573374595371563  sparsity=0.040388548057259666
```
```text
B_local=0.01  L2=0.13232747836123407  sparsity=0.013860486253124304
```
```text
B_local=0.0  L2=0.13212429667373482  sparsity=0.0
```

Required:

```text
nonincreasing error as B_local is tightened
strict measurable improvement from coarse to zero budget
```

Both pass.

## Growth-trigger sensitivity

```text
threshold=0.05  final_basis=10  enrich=[]  L2=0.19119002020634157
```
```text
threshold=0.035  final_basis=12  enrich=[70, 120]  L2=0.18008929313945907
```
```text
threshold=0.03  final_basis=13  enrich=[10, 20, 70]  L2=0.13820391719791275
```
```text
threshold=0.025  final_basis=13  enrich=[10, 20, 30]  L2=0.13232747836123407
```
```text
threshold=0.015  final_basis=13  enrich=[10, 20, 30]  L2=0.13232747836123407
```

The error falls by

```text
30.79 %
```

from the most restrictive to the converged trigger regime.

## Sparse audit requirements

```text
sampled audit failures         = 0
dense sentinel failures        = 0
dense sentinels                = 2
candidate peak-memory reduction >= 95%
dense audit pair-work reduction vs v0.17 >= 65%
```

Measured:

```text
sampled failures:
0

dense sentinels:
2

candidate memory reduction:
97.55 %

dense audit pair reduction:
71.15 %
```

## Release status

```text
passed = True
```

Every configured v0.18 acceptance check passes.
