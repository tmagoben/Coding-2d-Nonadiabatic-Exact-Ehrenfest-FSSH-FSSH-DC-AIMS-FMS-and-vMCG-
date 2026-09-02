# v0.21 Benchmark Campaign

## Scientific question

v0.21 asks whether the Gaussian dynamics core can be made **complex,
representation-neutral, arbitrary-state, and block sparse** before any physical spin
Hamiltonian is introduced.

The campaign deliberately avoids SOC. It stresses the mathematical structures that an
SOC backend would later require while keeping the framework equally valid for ordinary
spin-free nonadiabatic dynamics.

## Campaign layers

### 1. Complex electronic-operator covariance

A spin-free molecular LVC provider is lifted into the generic operator contract

$$
H(q),\qquad
K_a(q)=\langle\Phi|\partial_a\hat H_e|\Phi\rangle,\qquad
D_a(q)=\langle\Phi|\partial_a\Phi\rangle.
$$

A smooth coordinate-dependent complex unitary transformation

$$
G(q)\in U(2)
$$

is applied analytically.

The campaign verifies

$$
H'=G^\dagger HG,
$$

$$
K_a'=G^\dagger K_aG,
$$

and

$$
D_a'=G^\dagger D_aG+G^\dagger\partial_aG.
$$

Measured:

```text
H relative error:
0.0

maximum K_a relative error:
0.0

force error:
8.673617379884035e-19
```

### 2. Full Gaussian-block covariance

Three moving nuclear Gaussians carry complete 2-state electronic blocks.

The dense block matrices are constructed in two equivalent complex electronic gauges.

The expected transformation is

$$
S'=\mathcal G^\dagger S\mathcal G,
$$

$$
H'=\mathcal G^\dagger H\mathcal G,
$$

and

$$
T'
=
\mathcal G^\dagger T\mathcal G
+
\mathcal G^\dagger S\dot{\mathcal G}.
$$

Measured:

```text
S relative error:
1.7399964775418827e-16

H relative error:
1.7933278321452532e-16

T relative error:
1.7899468791987713e-16

maximum sparse-edge score change:
5.551115123125783e-17
```

### 3. Time-dependent gauge-equivalent propagation

The same prescribed moving Gaussian basis is propagated in the base electronic frame
and in the smoothly varying complex gauge.

The transformed initial coefficients satisfy

$$
C'_0=\mathcal G_0^\dagger C_0.
$$

At the final time the transformed solution is mapped back using

$$
C_{mapped}(t_f)=\mathcal G(t_f)C'(t_f).
$$

The comparison uses the physical final overlap metric.

| dt | Steps | Gauge-mapped coefficient error | Maximum norm drift |
|---:|---:|---:|---:|
| 0.02000 | 5 | 2.04745034e-08 | 6.57036980e-09 |
| 0.01000 | 10 | 5.11867320e-09 | 1.64260239e-09 |
| 0.00500 | 20 | 1.27967170e-09 | 4.10651957e-10 |

Observed orders:

```text
[1.999986653128339, 1.9999961654805536]
```

This campaign therefore checks the complete time-dependent inhomogeneous gauge
transformation, not only a static matrix similarity transform.

### 4. Degenerate-subspace stress test

An 8-state manifold is rotated by an arbitrary unitary matrix.

Instead of assigning individual roots, the entire overlap matrix is aligned using its
SVD/Procrustes solution.

```text
minimum singular value:
0.9999999999999998

anti-Hermitian residual:
1.967477117880543e-15
```

### 5. Wilson-loop gauge test

A 4-state complex electronic cycle is transformed by independent local unitary gauges.

The Wilson-loop eigenphases are required to remain invariant.

```text
maximum eigenphase change:
9.43689570931383e-16
```

### 6. Block sparse convergence

A 16-Gaussian, 4-state complex synthetic operator benchmark is compared with the
complete dense block reference.

Both the edge threshold and accumulated local omission budget are independently
relaxed.

The campaign requires non-increasing errors in all of

$$
S,\quad H,\quad T.
$$

At zero local omission budget:

```text
S error:
1.3757894522835919e-07

H error:
1.5178432514327997e-06

T error:
1.0997753094998538e-06
```

### 7. Dynamic topology

Six Gaussian blocks move through one another and then separate.

The sparse graph must demonstrate both entry and deletion events.

```text
entered edges:
15

exited edges:
9

maximum active edges:
15

final active edges:
6
```

### 8. Arbitrary electronic dimension

The same block-sparse architecture is exercised at

$$
s=2,4,8.
$$

The Gaussian graph remains identical while scalar block storage scales as \(s^2\).

```text
H_nnz / s^2:
[114.0, 114.0, 114.0]

relative scaling error:
0.0
```

## Reproduction

Run:

```bash
python examples/110_recompute_v021_campaign.py
```

Canonical machine-readable output:

```text
results/v021_complex_block_framework_campaign.json
```

## Acceptance philosophy

v0.21 is accepted only when every configured covariance, convergence, topology,
subspace, Wilson-loop, norm, and state-dimension criterion passes.

No SOC-specific observable is part of the v0.21 acceptance contract.
