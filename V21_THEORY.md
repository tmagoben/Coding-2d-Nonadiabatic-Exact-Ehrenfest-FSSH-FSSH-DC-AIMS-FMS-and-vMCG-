# v0.21 Theory

v0.21 generalizes the molecular sparse framework to a **complex, arbitrary-state,
block-valued electronic representation** without introducing spin physics.

The electronic contract is

$$
H(q)=H^\dagger(q),
$$

$$
K_a(q)=K_a^\dagger(q),
$$

and

$$
D_a(q)=-D_a^\dagger(q),
$$

where \(K_a\) contains matrix elements of the derivative of the physical electronic
operator.

Under \(G(q)\in U(s)\),

$$
H'=G^\dagger HG,
$$

$$
K_a'=G^\dagger K_aG,
$$

$$
D_a'=G^\dagger D_aG+G^\dagger\partial_aG.
$$

The Gaussian basis is block-complete locally:

$$
|\Psi_G\rangle
=
\sum_{i=1}^N
\sum_{\alpha=1}^s
C_{i\alpha}
|g_i\rangle|\phi_{i\alpha}\rangle.
$$

The coefficient dimension is \(Ns\).

Release covariance errors are:

```text
point H:
0.0

point dH:
0.0

block S:
1.7399964775418827e-16

block H:
1.7933278321452532e-16

block T:
1.7899468791987713e-16
```

The v0.21 sparse edge score is constructed only from gauge-invariant block quantities.
The maximum score change under the release complex gauge is

```text
5.551115123125783e-17
```

Detailed theory is in:

- `docs/01_MATHEMATICAL_FOUNDATIONS.md`
- `docs/02_COMPLEX_GAUGE_FRAMEWORK.md`
- `docs/03_BLOCK_SPARSE_ALGORITHM.md`
