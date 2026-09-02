# Theoretical Background

This file is retained for compatibility with the earlier repository layout.

The canonical derivation has been rewritten for v0.21 and now continues below. Related
documents are in `docs/`, especially `02_COMPLEX_GAUGE_FRAMEWORK.md` and
`03_BLOCK_SPARSE_ALGORITHM.md`.

This document derives the core equations used by the Gaussian dynamics framework.

## 1. Molecular time-dependent Schrödinger equation

Let electronic coordinates be $r$ and nuclear coordinates be $R$. In atomic units,

$$
i\frac{\partial}{\partial t}\Psi(r,R,t)
=
\hat H\Psi(r,R,t),
$$

with

$$
\hat H
=
\hat T_N+\hat H_e(R).
$$

The nuclear kinetic operator is

$$
\hat T_N
=
-\frac12
\sum_A
\frac{1}{M_A}
\nabla_A^2.
$$

No adiabatic, diabatic, spin, or SOC representation has yet been chosen.

## 2. Coordinate-dependent electronic basis

Choose $s$ electronic functions collected as columns of

$$
\Phi(R)
=
\begin{bmatrix}
|\phi_1(R)\rangle & \cdots & |\phi_s(R)\rangle
\end{bmatrix}.
$$

Assume

$$
\Phi^\dagger\Phi=I_s.
$$

Expand the molecular state as

$$
|\Psi\rangle
=
\sum_{\alpha=1}^s
\chi_\alpha(R,t)
|\phi_\alpha(R)\rangle.
$$

Define the matrix-valued derivative connection

$$
\boxed{
D_A(R)
=
\Phi^\dagger(R)
\nabla_A\Phi(R).
}
$$

Differentiating $\Phi^\dagger\Phi=I$ gives

$$
(\nabla_A\Phi^\dagger)\Phi
+
\Phi^\dagger\nabla_A\Phi
=
0,
$$

so

$$
\boxed{
D_A^\dagger=-D_A.
}
$$

The connection is anti-Hermitian.

## 3. Covariant nuclear derivative

For the electronic coefficient vector $\chi$,

$$
\nabla_A(\Phi\chi)
=
(\nabla_A\Phi)\chi
+
\Phi\nabla_A\chi.
$$

Projecting with $\Phi^\dagger$,

$$
\Phi^\dagger\nabla_A(\Phi\chi)
=
D_A\chi+\nabla_A\chi.
$$

Thus the natural derivative is

$$
\boxed{
\nabla_A^{\mathrm{cov}}
=
\nabla_A+D_A.
}
$$

For a complete modeled electronic subspace, the projected kinetic operator can be
written schematically as

$$
\boxed{
\hat T_N^{\mathrm{proj}}
=
-\frac12
\sum_A
\frac1{M_A}
(\nabla_A+D_A)^2.
}
$$

## 4. Local electronic gauge freedom

Any local unitary transformation

$$
G(R)\in U(s)
$$

defines an equally valid electronic basis

$$
\Phi'(R)=\Phi(R)G(R).
$$

The coefficient vector transforms oppositely,

$$
\boxed{
\chi'=G^\dagger\chi,
}
$$

so the physical wavefunction is unchanged.

The electronic Hamiltonian matrix transforms as

$$
\boxed{
H_e'=G^\dagger H_eG.
}
$$

The derivative connection transforms as

$$
\begin{aligned}
D_A'
&=
(\Phi G)^\dagger
\nabla_A(\Phi G)
\\
&=
G^\dagger D_AG
+
G^\dagger\nabla_AG.
\end{aligned}
$$

Therefore

$$
\boxed{
D_A'
=
G^\dagger D_AG
+
G^\dagger\nabla_AG.
}
$$

This is the non-Abelian connection transformation law.

## 5. Physical Hamiltonian-derivative operator

The framework distinguishes:

1. the derivative of the **matrix representation** $H_e(R)$;
2. the matrix elements of the derivative of the **physical operator**
   $\partial_A\hat H_e$.

v0.21 stores the second object:

$$
\boxed{
K_A
=
\Phi^\dagger
(\partial_A\hat H_e)
\Phi.
}
$$

It transforms homogeneously:

$$
\boxed{
K_A'=G^\dagger K_AG.
}
$$

In a nondegenerate adiabatic basis,

$$
H_e=\operatorname{diag}(E_1,\ldots,E_s).
$$

For $i\ne j$,

$$
\boxed{
(K_A)_{ij}
=
(E_j-E_i)(D_A)_{ij}.
}
$$

The diagonal elements are

$$
(K_A)_{ii}=\partial_AE_i.
$$

This relation follows by differentiating the electronic eigenvalue equation and
projecting between different eigenstates.

## 6. Gauge-invariant force expectation

For a normalized local electronic vector $c$,

$$
c^\dagger c=1,
$$

define

$$
\boxed{
F_A
=
-c^\dagger K_Ac.
}
$$

Under

$$
c'=G^\dagger c,
\qquad
K_A'=G^\dagger K_AG,
$$

the force is unchanged:

$$
F_A'=F_A.
$$

This gives a representation-neutral route to future mean-field or generalized nuclear
guidance policies.

## 7. Moving Gaussian nuclear basis

Use nuclear Gaussians

$$
g_i(R,t)
=
\mathcal N_i
\exp\left[
-\frac12(R-q_i)^TA_i(R-q_i)
+i\,p_i^T(R-q_i)
\right].
$$

v0.21 associates each nuclear Gaussian with a full local electronic block:

$$
\boxed{
|\Xi_{i\alpha}(t)\rangle
=
|g_i(t)\rangle
|\phi_{i\alpha}(t)\rangle.
}
$$

The approximate molecular state is

$$
\boxed{
|\Psi_G(t)\rangle
=
\sum_{i=1}^{N}
\sum_{\alpha=1}^{s}
C_{i\alpha}(t)
|\Xi_{i\alpha}(t)\rangle.
}
$$

The coefficient dimension is

$$
M=Ns.
$$

## 8. Projected moving-basis TDSE

Define

$$
S_{\mu\nu}
=
\langle\Xi_\mu|\Xi_\nu\rangle,
$$

$$
H_{\mu\nu}
=
\langle\Xi_\mu|\hat H|\Xi_\nu\rangle,
$$

and

$$
T_{\mu\nu}
=
\langle\Xi_\mu|\dot\Xi_\nu\rangle.
$$

Insert the basis expansion into the TDSE:

$$
i\frac{d}{dt}
\sum_\nu
C_\nu|\Xi_\nu\rangle
=
\hat H
\sum_\nu
C_\nu|\Xi_\nu\rangle.
$$

Expanding the time derivative,

$$
i\sum_\nu
\dot C_\nu|\Xi_\nu\rangle
+
i\sum_\nu
C_\nu|\dot\Xi_\nu\rangle
=
\sum_\nu
C_\nu\hat H|\Xi_\nu\rangle.
$$

Projecting with $\langle\Xi_\mu|$ gives

$$
iS\dot C+iTC=HC.
$$

Hence

$$
\boxed{
iS\dot C=(H-iT)C.
}
$$

Equivalently,

$$
\boxed{
S\dot C=-(iH+T)C.
}
$$

## 9. Generalized norm and metric compatibility

The represented-state norm is

$$
\boxed{
\langle\Psi_G|\Psi_G\rangle
=
C^\dagger SC.
}
$$

It is not generally $C^\dagger C$.

Differentiating $S_{\mu\nu}=\langle\Xi_\mu|\Xi_\nu\rangle$ gives

$$
\boxed{
\dot S=T+T^\dagger.
}
$$

This identity is the central metric-compatibility condition for a moving nonorthogonal
basis.

## 10. Midpoint/Cayley propagation

Let

$$
K=iH+T.
$$

A midpoint discretization of

$$
S\dot C=-KC
$$

gives

$$
S_m
\frac{C_{n+1}-C_n}{\Delta t}
=
-\frac12K_m(C_{n+1}+C_n).
$$

Rearranging,

$$
\boxed{
\left[
S_m+\frac{\Delta t}{2}K_m
\right]C_{n+1}
=
\left[
S_m-\frac{\Delta t}{2}K_m
\right]C_n.
}
$$

The implementation solves this linear system directly.

It never forms $S^{-1}$.

## 11. Why this framework can support SOC without requiring SOC

All equations above are valid for:

- real adiabatic electronic states;
- complex quasi-diabatic states;
- arbitrary smooth $U(s)$ electronic frames;
- a future Hamiltonian containing spin-orbit coupling.

A later SOC backend can supply a different complex Hermitian $H$ while leaving the
Gaussian propagation algebra unchanged.

That is the intended abstraction boundary.
