# Validation strategy

The project uses layered validation so that a successful scalar observable cannot hide
a representation, propagation, or sparsity failure.

## 1. Algebraic invariants

Require:

$$
H=H^\dagger,
\qquad
K_a=K_a^\dagger,
\qquad
D_a=-D_a^\dagger.
$$

Then verify exact complex gauge transformation laws.

```text
H covariance:
0.0

dH covariance:
0.0

force covariance:
8.673617379884035e-19
```

## 2. Block matrix covariance

Verify

$$
S'=\mathcal G^\dagger S\mathcal G,
$$

$$
H'=\mathcal G^\dagger H\mathcal G,
$$

$$
T'
=
\mathcal G^\dagger T\mathcal G
+
\mathcal G^\dagger S\dot{\mathcal G}.
$$

```text
S error:
1.7399964775418827e-16

H error:
1.7933278321452532e-16

T error:
1.7899468791987713e-16
```

## 3. Sparse score representation invariance

Maximum edge-score change under the complex gauge:

```text
5.551115123125783e-17
```

## 4. Time-dependent gauge propagation

| dt | Gauge-mapped coefficient error | Norm drift |
|---:|---:|---:|
| 0.02000 | 2.047450e-08 | 6.570370e-09 |
| 0.01000 | 5.118673e-09 | 1.642602e-09 |
| 0.00500 | 1.279672e-09 | 4.106520e-10 |

Observed orders:

```text
[1.999986653128339, 1.9999961654805536]
```

## 5. Degenerate-subspace robustness

```text
minimum singular value:
0.9999999999999998

aligned-overlap anti-Hermitian residual:
1.967477117880543e-15
```

No individual root assignment is needed.

## 6. Wilson-loop spectrum

Maximum eigenphase change under local complex gauges:

```text
9.43689570931383e-16
```

## 7. Sparse convergence

Both edge threshold and omitted-score budget are swept against a dense block reference.

The release requires non-increasing $S,H,T$ errors as the sparse approximation is
relaxed.

## 8. Dynamic graph topology

The graph must enter and exit edges during the stress trajectory.

```text
entered:
15

exited:
9
```

## 9. Electronic dimension

The same machinery is tested at

$$
s=2,4,8.
$$

No v0.21 core path assumes exactly two states.

## 10. Release philosophy

A release passes only if every configured acceptance condition passes.

No claim of full-wavefunction or representation correctness is inferred merely from
population agreement.
