# v0.13 Detailed Derivations

## A. Orthogonal projection residual

Let the current spinor-complete approximation space be

$$
\mathcal V_N
=
\operatorname{span}
\left\{
g_i(R)|d_a\rangle
:
i=1,\ldots,N,\;
a=0,1
\right\}.
$$

Let

$$
P_N
$$

be the orthogonal projector onto this space.

The current best approximation is

$$
|\Psi_N\rangle
=
P_N|\Psi\rangle.
$$

Define

$$
|r_N\rangle
=
|\Psi\rangle-|\Psi_N\rangle.
$$

For every

$$
|v\rangle\in\mathcal V_N,
$$

orthogonal projection gives

$$
\boxed{
\langle v|r_N\rangle=0.
}
$$

By the Pythagorean theorem,

$$
\boxed{
\|\Psi\|^2
=
\|\Psi_N\|^2
+
\|r_N\|^2.
}
$$

For normalized target states,

$$
\|\Psi\|^2=1,
$$

so

$$
\|r_N\|^2
=
1-\|\Psi_N\|^2.
$$

For an exact projection,

$$
F_N
=
\frac{
|\langle\Psi|\Psi_N\rangle|^2
}{
\|\Psi\|^2\|\Psi_N\|^2
}
=
\frac{\|\Psi_N\|^2}{\|\Psi\|^2}.
$$

Therefore

$$
\boxed{
F_N
=
1-
\frac{\|r_N\|^2}{\|\Psi\|^2}.
}
$$

---

## B. Nuclear candidate orthogonalization

Let

$$
\mathcal B_N
=
\operatorname{span}
\{g_1,\ldots,g_N\}.
$$

For normalized candidate

$$
g_c,
$$

write its projection as

$$
P_{\mathcal B_N}g_c
=
\sum_i\alpha_i g_i.
$$

The least-squares condition is

$$
\langle g_j|
g_c-\sum_i\alpha_i g_i
\rangle
=
0.
$$

Hence

$$
s_j
-
\sum_iS_{ji}\alpha_i
=
0,
$$

where

$$
S_{ji}
=
\langle g_j|g_i\rangle
$$

and

$$
s_j
=
\langle g_j|g_c\rangle.
$$

Therefore

$$
\boxed{
S\alpha=s.
}
$$

Assuming a nonsingular overlap matrix,

$$
\boxed{
\alpha=S^{-1}s.
}
$$

The orthogonalized candidate is

$$
\boxed{
g_c^\perp
=
g_c-\sum_i\alpha_i g_i.
}
$$

---

## C. Norm of the new direction

Compute

$$
\langle g_c^\perp|g_c^\perp\rangle.
$$

Expand:

$$
\begin{aligned}
\|g_c^\perp\|^2
&=
\langle g_c|g_c\rangle
-\alpha^\dagger s
-s^\dagger\alpha
+\alpha^\dagger S\alpha.
\end{aligned}
$$

Because

$$
S\alpha=s,
$$

$$
\alpha^\dagger S\alpha
=
\alpha^\dagger s.
$$

For Hermitian

$$
S,
$$

$$
\alpha^\dagger s
=
s^\dagger S^{-1}s
$$

is real and nonnegative.

With normalized

$$
g_c,
$$

$$
\boxed{
\|g_c^\perp\|^2
=
1-s^\dagger S^{-1}s.
}
$$

Define

$$
\boxed{
n_c
=
1-s^\dagger S^{-1}s.
}
$$

If

$$
n_c=0,
$$

the candidate lies in the existing span.

---

## D. Spinor-complete candidate subspace

The electronic basis is orthonormal:

$$
\langle d_a|d_b\rangle
=
\delta_{ab}.
$$

The new independent subspace introduced by one nuclear candidate is

$$
\mathcal W_c
=
\operatorname{span}
\left\{
g_c^\perp|d_0\rangle,
g_c^\perp|d_1\rangle
\right\}.
$$

The two basis vectors are orthogonal to one another electronically and both have norm

$$
n_c.
$$

---

## E. Optimal one-step residual coefficient

Write the electronic residual components

$$
|r\rangle
=
r_0(R)|d_0\rangle
+
r_1(R)|d_1\rangle.
$$

Approximate the residual component inside

$$
\mathcal W_c
$$

as

$$
|\delta\Psi_c\rangle
=
\sum_a
\beta_a
g_c^\perp|d_a\rangle.
$$

Minimize

$$
\left\|
r-\delta\Psi_c
\right\|^2
$$

with respect to

$$
\beta_a^*.
$$

For each electronic component,

$$
\frac{\partial}{\partial\beta_a^*}
\left\|
r_a-\beta_ag_c^\perp
\right\|^2
=
-
\langle g_c^\perp|r_a\rangle
+
\beta_a n_c.
$$

Set to zero:

$$
\boxed{
\beta_a
=
\frac{
\langle g_c^\perp|r_a\rangle
}{
n_c
}.
}
$$

---

## F. Why the old-basis projection term disappears

Because the current residual is orthogonal to every existing spinor-complete basis
function,

$$
\langle g_i d_a|r\rangle=0.
$$

Therefore

$$
\begin{aligned}
\langle g_c^\perp|r_a\rangle
&=
\left\langle
g_c-\sum_i\alpha_i g_i
\middle|
r_a
\right\rangle
\\
&=
\langle g_c|r_a\rangle
-
\sum_i\alpha_i^*
\langle g_i|r_a\rangle
\\
&=
\boxed{
\langle g_c|r_a\rangle.
}
\end{aligned}
$$

This is why the residual score is inexpensive once the current projection has been
computed.

---

## G. Exact squared residual reduction

For one component,

$$
\left\|
\beta_ag_c^\perp
\right\|^2
=
|\beta_a|^2n_c.
$$

Substitute

$$
\beta_a
=
\frac{\langle g_c^\perp|r_a\rangle}{n_c}.
$$

Then

$$
\left\|
\beta_ag_c^\perp
\right\|^2
=
\frac{
|\langle g_c^\perp|r_a\rangle|^2
}{
n_c
}.
$$

Sum over electronic components:

$$
\boxed{
\Delta_c
=
\frac{
\displaystyle\sum_a
|\langle g_c|r_a\rangle|^2
}{
n_c
}.
}
$$

Because the new correction is an orthogonal projection of the residual,

$$
\boxed{
\|r_{N+1}\|^2
=
\|r_N\|^2
-
\Delta_c.
}
$$

This is not an empirical score.

It is the exact one-candidate projection gain for the current residual and candidate
subspace.

---

## H. Relative residual reduction

The repository reports

$$
\epsilon_N
=
\frac{\|r_N\|^2}{\|\Psi\|^2}.
$$

Therefore the predicted decrease in the reported relative residual is

$$
\boxed{
\Delta\epsilon_c
=
\frac{\Delta_c}{\|\Psi\|^2}.
}
$$

The v0.13 regression suite compares this value with the actual reprojection result.

---

## I. Finite-grid vectorized form

Let the flattened candidate Gaussian values be a row vector

$$
G_c.
$$

Let the current residual be the matrix

$$
R
=
\begin{pmatrix}
r_0 & r_1
\end{pmatrix},
$$

with grid points along the row dimension.

Then

$$
\boxed{
b_c
=
G_c^\dagger R\,\Delta A
}
$$

is a two-component electronic residual-overlap vector.

For all candidates simultaneously, stack the Gaussian rows into

$$
G_{\mathrm{dict}}.
$$

Then

$$
\boxed{
B
=
G_{\mathrm{dict}}^\dagger R\,\Delta A
}
$$

is evaluated by one dense matrix multiplication.

This is the vectorized implementation used by the release benchmark.

---

## J. Candidate orthogonal norms for all dictionary entries

Let

$$
G_B
$$

contain the current nuclear basis functions.

The finite-grid nuclear overlap is

$$
S=G_B^\dagger G_B\,\Delta A.
$$

Let

$$
X
=
G_B^\dagger G_{\mathrm{dict}}\,\Delta A.
$$

Each candidate column

$$
x_c
$$

contains

$$
\langle g_i|g_c\rangle.
$$

Then

$$
\boxed{
n_c
=
\|g_c\|^2
-
x_c^\dagger S^+x_c.
}
$$

The implementation uses a pseudoinverse only for numerical robustness.

Candidates with

$$
n_c
$$

below the configured floor are rejected.

---

## K. Observable-aware second-stage screening

Residual gain determines the admissible short list.

Suppose

$$
\mathcal C_K
$$

contains the top

$$
K
$$

candidates by

$$
\Delta_c.
$$

For each candidate, the one-step corrected wavefunction is

$$
\boxed{
\Psi_c
=
\Psi_N
+
g_c^\perp
\begin{pmatrix}
\beta_0\\
\beta_1
\end{pmatrix}.
}
$$

Its reduced electronic density is

$$
\rho_c
=
\operatorname{Tr}_N
|\Psi_c\rangle\langle\Psi_c|.
$$

The benchmark selects

$$
\boxed{
c^*
=
\arg\min_{c\in\mathcal C_K}
\|\rho_c-\rho_{\mathrm{target}}\|_F.
}
$$

This uses only the known initial target state.

---

## L. Why the screening is not future fitting

The selection functional depends on

$$
\rho_{\mathrm{target}}(0),
$$

not

$$
\rho_{\mathrm{target}}(t_f).
$$

The exact final-time benchmark is not evaluated until after the complete initial bank
has been selected.

Therefore there is no leakage of future exact dynamics into basis construction.

---

## M. Time-dependent Galerkin defect

For the moving approximation

$$
|\Psi_G(t)\rangle
=
\sum_\mu
C_\mu(t)
|\Xi_\mu(t)\rangle,
$$

the exact TDSE is

$$
i|\dot\Psi\rangle
=
\hat H|\Psi\rangle.
$$

Define

$$
\boxed{
|\mathcal R\rangle
=
i|\dot\Psi_G\rangle
-
\hat H|\Psi_G\rangle.
}
$$

The defect norm is

$$
\boxed{
\eta
=
\|\mathcal R\|.
}
$$

A normalized diagnostic used by the code is

$$
\boxed{
\eta_H
=
\frac{\|\mathcal R\|}{\|\hat H\Psi_G\|}.
}
$$

---

## N. Coefficient derivative

The moving-basis equation is

$$
iS\dot C
=
(H-iT)C.
$$

Rearrange:

$$
iS\dot C
=
HC-iTC.
$$

Multiply by

$$
-i:
$$

$$
S\dot C
=
-iHC-TC.
$$

Therefore

$$
\boxed{
\dot C
=
S^{-1}
[-(iH+T)C].
}
$$

The implementation solves the linear system rather than forming

$$
S^{-1}
$$

explicitly.

---

## O. Total wavefunction derivative

For spinor-complete global-diabatic Gaussians,

$$
\Psi_G(R,t)
=
\sum_i
g_i(R,t)\mathbf C_i(t).
$$

Therefore

$$
\boxed{
\dot\Psi_G
=
\sum_i
\dot g_i\mathbf C_i
+
\sum_i
g_i\dot{\mathbf C}_i.
}
$$

Both terms are required.

Dropping the moving-Gaussian term would produce an incorrect TDSE residual.

---

## P. Applying the exact benchmark Hamiltonian

The diagnostic Hamiltonian is

$$
\hat H
=
-\frac{1}{2M}\nabla^2I
+
V_d(R).
$$

The kinetic term is applied on the same periodic FFT grid as the exact benchmark:

$$
\widetilde{T\Psi}(k)
=
\frac{k_x^2+k_y^2}{2M}
\widetilde{\Psi}(k).
$$

The potential action is

$$
\boxed{
[V\Psi](R)
=
V_d(R)\Psi(R).
}
$$

This keeps the TDSE-defect diagnostic independent of the finite Gaussian matrix
projection.

---

## Q. Galerkin orthogonality

Let the projected equations be solved exactly in the current basis.

Then for every current tangent/basis test direction

$$
|\delta\Psi\rangle,
$$

the Dirac-Frenkel/Galerkin condition is

$$
\boxed{
\langle\delta\Psi|\mathcal R\rangle=0
}
$$

for the directions included in the projected equations.

v0.13 numerically projects

$$
\mathcal R
$$

back onto the current spinor-complete basis.

The projected norm is required to be much smaller than the total defect norm.

---

## R. Defect-driven candidate gain

Unlike the initial projection residual, the numerically evaluated TDSE defect need not
be exactly orthogonal to the candidate nuclear span on the finite diagnostic grid.

Therefore v0.13 explicitly orthogonalizes the candidate:

$$
g_c^\perp
=
g_c-\sum_i\alpha_i g_i.
$$

Define

$$
b_a
=
\langle g_c^\perp|\mathcal R_a\rangle.
$$

The squared defect captured by the candidate pair is

$$
\boxed{
\Delta_c^{\mathrm{TDSE}}
=
\frac{
\sum_a|b_a|^2
}{
n_c
}.
}
$$

---

## S. Zero-coefficient enrichment and wavefunction continuity

Let the old coefficient vector be

$$
C.
$$

After adding a complete two-component electronic pair, define

$$
\boxed{
C'
=
\begin{pmatrix}
C\\
0\\
0
\end{pmatrix}.
}
$$

Then

$$
\boxed{
\Psi'(R)=\Psi(R).
}
$$

However, the new coefficient derivative

$$
\dot C'
$$

is solved in the larger Galerkin space.

Thus

$$
\mathcal R'
$$

can be smaller even though the instantaneous wavefunction is unchanged.

---

## T. Predicted defect reduction

If the newly added pair is the only new orthogonal tangent direction and the Galerkin
solve is exact, then

$$
\boxed{
\|\mathcal R'\|^2
=
\|\mathcal R\|^2
-
\Delta_c^{\mathrm{TDSE}}.
}
$$

The release benchmark obtains

$$
\Delta_c^{\mathrm{predicted}}
=
0.0171436285450576,
$$

and

$$
\Delta_c^{\mathrm{actual}}
=
0.017143629785279.
$$

Their relative difference is at the level of floating-point/grid numerical error.

---

## U. Why a large residual gain is not enough by itself

Suppose a candidate has

$$
n_c\ll1.
$$

Then it is almost contained in the current span.

The gain expression contains

$$
1/n_c,
$$

so numerical roundoff can become strongly amplified.

Therefore the algorithm also requires:

1. a minimum orthogonal norm;
2. a maximum expanded overlap condition number.

Residual reduction and conditioning are treated as simultaneous numerical constraints.

---

## V. Residual monotonicity versus observable monotonicity

Pure Hilbert residual greedy guarantees

$$
\boxed{
\|r_{N+1}\|^2
\le
\|r_N\|^2.
}
$$

It does **not** imply that every nonlinear reduced observable improves monotonically.

For example,

$$
\|\rho_N-\rho_{\mathrm{target}}\|_F
$$

may temporarily increase even while the full wavefunction residual decreases.

This is why v0.13 reports both quantities separately and does not claim that one is a
proxy for the other.

---

## W. Error hierarchy after v0.13

For the release reference,

$$
\epsilon_{\mathrm{init}}
\approx
3.21\times10^{-2},
$$

while the representation-consistent projected dynamics error is only

$$
\epsilon_{\mathrm{dyn}}
\approx
1.14\times10^{-4}.
$$

Thus

$$
\boxed{
\epsilon_{\mathrm{dyn}}
\ll
\epsilon_{\mathrm{init}}.
}
$$

The dominant error remains finite representation of the intended initial state, but
v0.13 now reduces that error through a measurable and reproducible adaptive criterion
rather than hand-selected basis placement.
