# Complex gauge framework

v0.21 is representation neutral and deliberately does **not** introduce spin physics.

For an orthonormal electronic frame \(\Phi(q)\), define

$$
H=\Phi^\dagger\hat H_e\Phi,
\qquad
K_a=\Phi^\dagger(\partial_a\hat H_e)\Phi,
\qquad
D_a=\Phi^\dagger\partial_a\Phi.
$$

The required algebraic properties are

$$
H=H^\dagger,
\qquad
K_a=K_a^\dagger,
\qquad
D_a=-D_a^\dagger.
$$

## Local U(s) transformation

For

$$
\Phi'=\Phi G(q),\qquad G(q)\in U(s),
$$

the coefficient vector changes as \(c'=G^\dagger c\). Physical operator matrices transform homogeneously,

$$
H'=G^\dagger HG,
\qquad
K_a'=G^\dagger K_aG,
$$

whereas the derivative connection transforms inhomogeneously,

$$
\boxed{D_a'=G^\dagger D_aG+G^\dagger\partial_aG.}
$$

This distinction is fundamental. The derivative of a matrix representation is not the same object as the matrix elements of the derivative of the physical operator.

## Gauge-invariant force

For a normalized local electronic vector,

$$
F_a=-c^\dagger K_ac.
$$

Under the transformation above,

$$
F_a'=-c^\dagger GG^\dagger K_aGG^\dagger c=F_a.
$$

## Gaussian block covariance

If every Gaussian carries a complete local electronic block, define

$$
\mathcal G=\operatorname{diag}(G_1,\ldots,G_N).
$$

Then

$$
S'=\mathcal G^\dagger S\mathcal G,
\qquad
H'=\mathcal G^\dagger H\mathcal G.
$$

Because the basis itself moves,

$$
\boxed{
T'=\mathcal G^\dagger T\mathcal G+\mathcal G^\dagger S\dot{\mathcal G}.
}
$$

The extra term means that \(\|T_{ij}\|\) is not a representation-invariant sparsification criterion under a coordinate-dependent gauge.

## Degenerate subspaces

Inside a degenerate manifold, individual eigenvectors are not unique. For the subspace overlap

$$
O=U\Sigma V^\dagger,
$$

v0.21 uses the Procrustes transform

$$
\boxed{W=VU^\dagger}
$$

so that

$$
OW=U\Sigma U^\dagger
$$

is Hermitian positive semidefinite. This aligns the full subspace without claiming a unique root identity.

## Wilson loops

For electronic links around a closed cycle,

$$
W=L_{01}L_{12}\cdots L_{m0}.
$$

Under local gauges \(L_{ij}'=G_i^\dagger L_{ij}G_j\),

$$
W'=G_0^\dagger WG_0.
$$

Therefore the Wilson eigenvalues/eigenphases are gauge invariant. Local gauge smoothing must not be confused with removal of genuine global holonomy.

## SOC readiness without SOC dependence

A later backend may provide

$$
H=H_{\rm spin\mbox{-}free}+H_{\rm SOC},
$$

but the same interface also accepts \(H_{\rm SOC}=0\). The Gaussian engine therefore remains a general nonadiabatic framework rather than becoming an SOC-specific codebase.
