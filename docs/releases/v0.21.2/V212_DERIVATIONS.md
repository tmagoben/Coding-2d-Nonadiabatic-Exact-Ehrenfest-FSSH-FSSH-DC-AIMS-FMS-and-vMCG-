# v0.21.2 Detailed Derivations — Pre-SOC Integration Hardening

v0.21.2 introduces no spin Hamiltonian. It hardens the representation-neutral complex
block framework so the first SOC term can later be added as an electronic-operator
backend rather than forcing another redesign of the Gaussian engine.

## 1. Unequal-width Gaussian pair algebra

For normalized real-width frozen Gaussians

$$
g_i(q)=N_i\exp\left[-\frac12(q-q_i)^TA_i(q-q_i)+ip_i^T(q-q_i)\right],
$$

with independent positive-definite widths $A_i$ and $A_j$, define

$$
B=A_i+A_j,
$$

and

$$
\ell=A_iq_i+A_jq_j+i(p_j-p_i).
$$

The exact overlap is

$$
\boxed{
S_{ij}^{\mathrm{nuc}}
=
\frac{2^{d/2}(\det A_i\det A_j)^{1/4}}
{\sqrt{\det(A_i+A_j)}}
\exp\left[c+\frac12\ell^TB^{-1}\ell\right],
}
$$

where

$$
c=-\frac12q_i^TA_iq_i-\frac12q_j^TA_jq_j
+ip_i^Tq_i-ip_j^Tq_j.
$$

The complex cross centroid is

$$
\mu=B^{-1}\ell,
$$

and its covariance is

$$
\Sigma=B^{-1}.
$$

The exact unequal-width kinetic matrix element used by v0.21.2 is

$$
\boxed{
T_{ij}^{\mathrm{nuc}}
=
\frac12S_{ij}^{\mathrm{nuc}}
\left[
 u_i^TM^{-1}u_j
+
\operatorname{tr}(A_iM^{-1}A_j\Sigma)
\right],
}
$$

with

$$
u_i=-A_i(\mu-q_i)-ip_i,
\qquad
u_j=-A_j(\mu-q_j)+ip_j.
$$

For molecular electronic evaluation a real geometry is required, so the pair centroid is
the real envelope saddle point

$$
\boxed{
q_c=(A_i+A_j)^{-1}(A_iq_i+A_jq_j).
}
$$

For equal widths this reduces to $(q_i+q_j)/2$.

The full block pair algebra remains

$$
S_{ij}=S_{ij}^{\mathrm{nuc}}O_{ij},
$$

$$
H_{ij}=T_{ij}^{\mathrm{nuc}}O_{ij}+S_{ij}^{\mathrm{nuc}}H_{ij}^{e},
$$

and therefore does not require a common Gaussian width.

## 2. Self-consistent local electronic guidance

Each Gaussian carries a local electronic coefficient block $c_i\in\mathbb{C}^s$.
For non-negligible local amplitude define

$$
\tilde c_i=\frac{c_i}{\sqrt{c_i^\dagger c_i}}.
$$

Let

$$
K_a(q_i)=\langle\Phi|\partial_a\hat H_e|\Phi\rangle_{q_i}.
$$

The generalized force policy is

$$
\boxed{
F_{i,a}
=-\tilde c_i^\dagger K_a(q_i)\tilde c_i.
}
$$

Under a local gauge

$$
c_i'=G_i^\dagger c_i,
\qquad
K_a'=G_i^\dagger K_aG_i,
$$

so

$$
F_{i,a}'=F_{i,a}.
$$

Thus the nuclear guidance policy is representation neutral and does not depend on spin
labels.

The current release uses a velocity-Verlet predictor/corrector:

$$
p_{i,n+1/2}=p_{i,n}+\frac{\Delta t}{2}F_{i,n},
$$

$$
q_{i,n+1}=q_{i,n}+\Delta t\,M_i^{-1}p_{i,n+1/2},
$$

followed by block coefficient propagation on the moved basis, re-evaluation of
$F_{i,n+1}$, and

$$
\boxed{
p_{i,n+1}=p_{i,n+1/2}+\frac{\Delta t}{2}F_{i,n+1}.
}
$$

A small fixed-point corrector makes the endpoint momentum, block coefficients, and force
mutually consistent.

This is a **general nuclear guidance policy**, not a claim of the full variational AIMS
nuclear equations.

## 3. Zero-block Gaussian birth

A new Gaussian adds an $s$-component coefficient block.
To leave the represented state unchanged at the instant of birth,

$$
\boxed{
c_{\mathrm{new}}=0_s.}
$$

If the old coefficient vector is $C$, insertion is

$$
C'=(C_1,\ldots,C_k,0_s,C_{k+1},\ldots)^T.
$$

Therefore the molecular wavefunction before and immediately after insertion is exactly
the same.

## 4. Metric-projected block pruning

Partition the overlap matrix into retained and deleted scalar components:

$$
S=
\begin{pmatrix}
S_{rr} & S_{rd}\\
S_{dr} & S_{dd}
\end{pmatrix},
\qquad
C=
\begin{pmatrix}C_r\\C_d\end{pmatrix}.
$$

The orthogonal projection of the represented state into the retained span obeys

$$
S_{rr}C_r'=S_{rr}C_r+S_{rd}C_d,
$$

so

$$
\boxed{
C_r'=C_r+S_{rr}^{-1}S_{rd}C_d.
}
$$

The squared lost norm is the Schur-complement expression

$$
\boxed{
\epsilon_{\mathrm{prune}}^2
=
C_d^\dagger
\left(
S_{dd}-S_{dr}S_{rr}^{-1}S_{rd}
\right)
C_d.
}
$$

The implementation solves the linear systems; it does not explicitly invert $S_{rr}$.

## 5. Representation-neutral electronic observables

For a physical Hermitian electronic operator $\hat O$, the local matrix is

$$
O(q)=\Phi^\dagger(q)\hat O\Phi(q).
$$

At a Gaussian pair centroid the electronic frames are transported into one common frame,
producing

$$
O_{ij}^{e}=U_{ci}^\dagger O(q_c)U_{cj}.
$$

The molecular block is

$$
\boxed{
\mathcal O_{ij}=S_{ij}^{\mathrm{nuc}}O_{ij}^{e}.
}
$$

The observable expectation is

$$
\boxed{
\langle O\rangle
=
\frac{C^\dagger\mathcal O C}{C^\dagger SC}.
}
$$

This is the future extension point for spin operators, charge operators, state-character
projectors, or any other electronic observable. No spin observable is built into the
core.

## 6. Full-subspace continuity diagnostics

The block framework already tolerates arbitrary local electronic frames, so v0.21.2
does **not** force root-by-root phase smoothing.

For a cross-geometry subspace overlap

$$
O=U\Sigma V^\dagger,
$$

its singular values measure subspace continuity. The Procrustes map

$$
W=VU^\dagger
$$

is recorded as a diagnostic alignment, while the provider's raw local frame is returned
unchanged.

This avoids inserting an artificial geometry-dependent gauge rotation whose derivative
would otherwise have to be added to the connection.

## 7. Pre-SOC complex-dtype invariant

The generalized core treats

$$
H,\quad K_a,\quad D_a,\quad O_{ij},\quad C
$$

as complex-valued quantities.

Real casts are allowed only for quantities that are physically real in this framework:

$$
q,\ p,\ M,\ F,
$$

plus the explicit legacy spin-free adiabatic adapter whose source contract is real.

The v0.21.2 audit classifies the latter separately rather than silently accepting an
unknown complex-to-real conversion.
