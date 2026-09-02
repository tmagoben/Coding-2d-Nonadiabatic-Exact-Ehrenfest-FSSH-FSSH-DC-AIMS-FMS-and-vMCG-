# v0.12 Detailed Derivations

---

## A. Reduced electronic density of a coordinate-dependent adiabatic packet

Let the full molecular state be

$$
|\Psi\rangle
=
\int dR\,
g(R)
|R\rangle
|\Phi_a(R)\rangle.
$$

Choose a fixed global diabatic electronic basis

$$
\{|d_\mu\rangle\}.
$$

Expand

$$
|\Phi_a(R)\rangle
=
\sum_\mu
U_{\mu a}(R)|d_\mu\rangle.
$$

Then

$$
|\Psi\rangle
=
\sum_\mu
\int dR\,
g(R)U_{\mu a}(R)
|R\rangle|d_\mu\rangle.
$$

The full density operator is

$$
|\Psi\rangle\langle\Psi|
=
\sum_{\mu\nu}
\int dR\,dR'\,
g(R)g^*(R')
U_{\mu a}(R)
U_{\nu a}^*(R')
|R\rangle\langle R'|
\otimes
|d_\mu\rangle\langle d_\nu|.
$$

Trace over nuclei:

$$
\operatorname{Tr}_R
[
|R\rangle\langle R'|
]
=
\delta(R-R').
$$

Therefore

$$
\boxed{
\rho_{\mu\nu}^{(e)}
=
\int dR\,
|g(R)|^2
U_{\mu a}(R)
U_{\nu a}^*(R).
}
$$

If $U_{\mu a}$ changes across the packet, $\rho_e$ is generally mixed.

---

## B. Purity of the center-frozen approximation

The center-frozen state is

$$
|\Psi_{\rm cf}\rangle
=
\int dR\,
g(R)|R\rangle
|\Phi_a(q_0)\rangle.
$$

Because the electronic vector is independent of $R$,

$$
|\Psi_{\rm cf}\rangle
=
|\chi\rangle_{\rm nuc}
\otimes
|\Phi_a(q_0)\rangle_{\rm el}.
$$

The electronic reduced density is

$$
\rho_e
=
|\Phi_a(q_0)\rangle
\langle\Phi_a(q_0)|.
$$

Then

$$
\rho_e^2=\rho_e.
$$

Hence

$$
\boxed{
\operatorname{Tr}(\rho_e^2)=1.
}
$$

Thus a center-frozen one-spinor basis can never represent reduced electronic
mixedness generated purely by the coordinate dependence of one adiabatic state.

---

## C. Spinor-complete Gaussian normal equations

Basis functions are

$$
|B_{ia}\rangle
=
g_i(R)|d_a\rangle.
$$

Approximate

$$
|\Psi\rangle
\approx
\sum_{ia}
C_{ia}|B_{ia}\rangle.
$$

Define residual

$$
|r\rangle
=
|\Psi_{\rm target}\rangle
-
\sum_{ia}C_{ia}|B_{ia}\rangle.
$$

Minimize

$$
\mathcal L
=
\langle r|r\rangle.
$$

Differentiate with respect to $C_{jb}^*$:

$$
\frac{\partial\mathcal L}{\partial C_{jb}^*}
=
-
\langle B_{jb}|\Psi_{\rm target}\rangle
+
\sum_{ia}
\langle B_{jb}|B_{ia}\rangle
C_{ia}.
$$

Set to zero:

$$
\sum_{ia}
S_{jb,ia}C_{ia}
=
b_{jb}.
$$

Hence

$$
\boxed{
SC=b.
}
$$

Because

$$
\langle d_b|d_a\rangle=\delta_{ba},
$$

$$
\boxed{
S_{jb,ia}
=
S_{ji}^{\rm nuc}\delta_{ba}.
}
$$

---

## D. Projection norm identity

For the exact least-squares projection,

$$
SC=b.
$$

The overlap of target and projection is

$$
\langle\Psi_{\rm target}|\Psi_{\rm proj}\rangle
=
b^\dagger C.
$$

Using

$$
b=SC,
$$

$$
b^\dagger C
=
C^\dagger S C.
$$

But

$$
C^\dagger SC
=
\langle\Psi_{\rm proj}|\Psi_{\rm proj}\rangle.
$$

Thus for an exact orthogonal projection,

$$
\boxed{
\langle\Psi_{\rm target}|\Psi_{\rm proj}\rangle
=
\|\Psi_{\rm proj}\|^2.
}
$$

Therefore

$$
\|\Psi_{\rm target}-\Psi_{\rm proj}\|^2
=
\|\Psi_{\rm target}\|^2
-
\|\Psi_{\rm proj}\|^2.
$$

In the finite grid implementation small differences can arise from finite-domain and
least-squares tolerances.

---

## E. Exact unequal-width Gaussian moments

For

$$
g_i^*g_j
\propto
\exp
\left[
-\frac12(q-\mu)^TB(q-\mu)
\right],
$$

with

$$
B=A_i+A_j,
$$

the covariance is

$$
\Sigma=B^{-1}.
$$

Therefore

$$
\boxed{
\langle q_\alpha\rangle_{ij}
=
\mu_\alpha S_{ij},
}
$$

and

$$
\boxed{
\langle q_\alpha q_\beta\rangle_{ij}
=
(
\mu_\alpha\mu_\beta
+
\Sigma_{\alpha\beta}
)
S_{ij}.
}
$$

These are algebraic complex moments of the cross density.

They are not moments of a positive probability distribution and do not use
$|\mu|^2$.

---

## F. Exact LVC potential integral

Write

$$
V_d
=
V_0I
+
V_x\sigma_x
+
V_z\sigma_z,
$$

with

$$
V_0
=
\frac12\omega^2(x^2+y^2),
$$

$$
V_x=\lambda y,
$$

$$
V_z=\kappa x.
$$

Then

$$
\langle g_i|V_d|g_j\rangle
=
\langle V_0\rangle I
+
\langle V_x\rangle\sigma_x
+
\langle V_z\rangle\sigma_z.
$$

Use the first and second moments:

$$
\langle V_0\rangle
=
\frac12\omega^2
S_{ij}
[
\mu_x^2+\mu_y^2+\Sigma_{xx}+\Sigma_{yy}
],
$$

$$
\langle V_x\rangle
=
\lambda\mu_yS_{ij},
$$

$$
\langle V_z\rangle
=
\kappa\mu_xS_{ij}.
$$

Therefore

$$
\boxed{
\begin{aligned}
V_{ij}^{\rm LVC}
=
S_{ij}
\Big[
&
\frac12\omega^2
(
\mu_x^2+\mu_y^2+
\Sigma_{xx}+\Sigma_{yy}
)I
\\
&+\lambda\mu_y\sigma_x
+\kappa\mu_x\sigma_z
\Big].
\end{aligned}
}
$$

This is exact for the model.

---

## G. Spinor-complete Hamiltonian block

Let

$$
T_{ij}^{\rm nuc}
=
\langle g_i|\hat T_N|g_j\rangle.
$$

For global diabatic basis states,

$$
\langle d_a|\hat T_N|d_b\rangle
=
\delta_{ab}\hat T_N.
$$

Thus

$$
\boxed{
H_{ia,jb}
=
T_{ij}^{\rm nuc}\delta_{ab}
+
[V_{ij}^{\rm LVC}]_{ab}.
}
$$

---

## H. Reduced density of the spinor-complete Gaussian wavefunction

The diabatic nuclear amplitudes are

$$
\psi_a(R)
=
\sum_i
C_{ia}g_i(R).
$$

Then

$$
\rho_{ab}
=
\int
\psi_a(R)\psi_b^*(R)dR.
$$

Expand:

$$
\rho_{ab}
=
\sum_{ij}
C_{ia}C_{jb}^*
\int
g_i(R)g_j^*(R)dR.
$$

Since

$$
\int g_i g_j^*
=
\langle g_j|g_i\rangle
=
S_{ji}^{\rm nuc},
$$

$$
\boxed{
\rho_{ab}
=
\sum_{ij}
C_{ia}C_{jb}^*
S_{ji}^{\rm nuc}.
}
$$

Matrix form:

$$
\boxed{
\rho
=
C^T
(S^{\rm nuc})^T
C^*.
}
$$

The implementation symmetrizes only accumulated floating-point anti-Hermiticity.

---

## I. Generalized norm of the spinor-complete basis

The total norm is

$$
\langle\Psi|\Psi\rangle
=
\sum_a
\sum_{ij}
C_{ia}^*
S_{ij}^{\rm nuc}
C_{ja}.
$$

Therefore

$$
\boxed{
N
=
\sum_a
C_a^\dagger
S^{\rm nuc}
C_a.
}
$$

This equals the trace of the unnormalized reduced electronic density.

---

## J. Moving-basis equation

For moving nonorthogonal basis functions,

$$
|\Psi\rangle
=
\sum_j
C_j(t)|\Xi_j(t)\rangle.
$$

Differentiate:

$$
|\dot\Psi\rangle
=
\sum_j
\dot C_j|\Xi_j\rangle
+
\sum_j
C_j|\dot\Xi_j\rangle.
$$

Insert into

$$
i|\dot\Psi\rangle
=
\hat H|\Psi\rangle.
$$

Project with $\langle\Xi_i|$:

$$
i
\sum_j
S_{ij}\dot C_j
+
i
\sum_j
T_{ij}C_j
=
\sum_j
H_{ij}C_j,
$$

where

$$
T_{ij}
=
\langle\Xi_i|\dot\Xi_j\rangle.
$$

Hence

$$
\boxed{
iS\dot C
=
(H-iT)C.
}
$$

---

## K. Metric identity

Differentiate

$$
S_{ij}
=
\langle\Xi_i|\Xi_j\rangle.
$$

Then

$$
\dot S_{ij}
=
\langle\dot\Xi_i|\Xi_j\rangle
+
\langle\Xi_i|\dot\Xi_j\rangle.
$$

Therefore

$$
\boxed{
\dot S
=
T^\dagger+T.
}
$$

The discrete v0.12 propagator enforces the endpoint analogue by minimally correcting
the Hermitian part of the midpoint seed.

---

## L. Midpoint Cayley equation

Write

$$
iS\dot C
=
(H-iT)C.
$$

Define

$$
K=iH+T.
$$

Then

$$
S\dot C=-KC.
$$

Midpoint discretization gives

$$
S_m
\frac{C_{n+1}-C_n}{\Delta t}
=
-
K_m
\frac{C_{n+1}+C_n}{2}.
$$

Rearrange:

$$
\boxed{
\left(
S_m+\frac{\Delta t}{2}K_m
\right)
C_{n+1}
=
\left(
S_m-\frac{\Delta t}{2}K_m
\right)
C_n.
}
$$

For

$$
S=I,
\qquad
T=0,
$$

this reduces to the standard Cayley/Crank-Nicolson propagator

$$
\boxed{
C_{n+1}
=
\left(I+\frac{i\Delta t}{2}H\right)^{-1}
\left(I-\frac{i\Delta t}{2}H\right)C_n.
}
$$

---

## M. Coordinate-dependent Born-Huang Gaussian basis

Define

$$
\Xi_i(R)
=
g_i(R)\Phi_{a_i}(R).
$$

On the two-dimensional benchmark grid, construct the diabatic vector field

$$
\boldsymbol\Xi_i(R)
=
g_i(R)
\mathbf U_{a_i}(R).
$$

The projected overlap is

$$
\boxed{
S_{ij}
=
\sum_R
\boldsymbol\Xi_i^\dagger(R)
\boldsymbol\Xi_j(R)
\Delta A.
}
$$

The projected Hamiltonian is

$$
\boxed{
H_{ij}
=
\sum_R
\boldsymbol\Xi_i^\dagger(R)
[
\hat T_{\rm FFT}
+
V_d(R)
]
\boldsymbol\Xi_j(R)
\Delta A.
}
$$

This avoids writing explicit second derivative couplings.

---

## N. Spectral kinetic action

For each electronic component,

$$
\tilde\Xi(k_x,k_y)
=
\mathcal F[\Xi(x,y)].
$$

The kinetic energy in momentum space is

$$
\boxed{
T(k_x,k_y)
=
\frac{k_x^2+k_y^2}{2M}.
}
$$

Therefore

$$
\boxed{
\hat T_{\rm FFT}\Xi
=
\mathcal F^{-1}
\left[
\frac{k_x^2+k_y^2}{2M}
\tilde\Xi
\right].
}
$$

This is the same momentum lattice used by the exact split-operator benchmark.

---

## O. Projected dynamics error

Let

$$
\rho_{\rm G}(t)
$$

be the Gaussian result initialized from the projected state.

Let

$$
\rho_{\rm proj}^{\rm exact}(t)
$$

be exact grid propagation from the identical projected initial wavefunction.

Define

$$
\boxed{
\epsilon_{\rm dyn}
=
\|
\rho_G(t)
-
\rho_{\rm proj}^{\rm exact}(t)
\|_F.
}
$$

For the nine-Gaussian v0.12 reference,

$$
\epsilon_{\rm dyn}
\approx
2.90\times10^{-4}.
$$

---

## P. Initial representation error

Let

$$
\rho_{\rm target}(0)
$$

be the intended exact initial density and

$$
\rho_{\rm proj}(0)
$$

the projected initial density.

Define

$$
\boxed{
\epsilon_{\rm init}
=
\|
\rho_{\rm proj}(0)
-
\rho_{\rm target}(0)
\|_F.
}
$$

For the nine-Gaussian reference,

$$
\epsilon_{\rm init}
\approx
3.545\times10^{-2}.
$$

---

## Q. Final target error

The final target error is

$$
\boxed{
\epsilon_{\rm target}
=
\|
\rho_G(t_f)
-
\rho_{\rm target}^{\rm exact}(t_f)
\|_F.
}
$$

For the nine-Gaussian reference,

$$
\epsilon_{\rm target}
\approx
3.500\times10^{-2}.
$$

The near equality

$$
\epsilon_{\rm target}
\approx
\epsilon_{\rm init}
$$

combined with

$$
\epsilon_{\rm dyn}\ll\epsilon_{\rm init}
$$

is the quantitative reason v0.12 identifies initial representation as the current
limiting approximation.

---

## R. Trace distance

For

$$
\Delta=\rho-\sigma,
$$

Hermitian with eigenvalues

$$
\lambda_k,
$$

the trace norm is

$$
\|\Delta\|_1
=
\sum_k|\lambda_k|.
$$

Therefore

$$
\boxed{
D_{\rm tr}
=
\frac12
\sum_k|\lambda_k|.
}
$$

The nine-Gaussian v0.12 reference gives approximately

$$
D_{\rm tr}
\approx
2.475\times10^{-2}.
$$

---

## S. Bloch vector

For

$$
\rho
=
\begin{pmatrix}
\rho_{00}&c\\
c^*&\rho_{11}
\end{pmatrix},
$$

write

$$
\rho
=
\frac12
(
I+r_x\sigma_x+r_y\sigma_y+r_z\sigma_z
).
$$

Matching matrix elements gives

$$
\boxed{
r_x=2\operatorname{Re}c,
}
$$

$$
\boxed{
r_y=-2\operatorname{Im}c,
}
$$

$$
\boxed{
r_z=\rho_{00}-\rho_{11}.
}
$$

Hence the Euclidean Bloch-vector error combines diagonal population imbalance and
off-diagonal coherence in one three-component diagnostic.

---

## T. Why adding spawned functions can fail to improve a projected initial bank

Suppose the initial nine-Gaussian bank already spans the relevant short-time nuclear
subspace well.

A new spawned Gaussian that is strongly overlapping with the existing span can:

1. increase $\kappa(S)$;
2. add little new projected Hilbert-space direction;
3. amplify coefficient sensitivity;
4. leave the physically important reduced-density error unchanged.

Therefore a spawn threshold alone is not a convergence criterion.

The relevant question is whether the new TBF reduces a measured residual or observable
error.

This motivates residual-driven basis enrichment as a natural later extension.
