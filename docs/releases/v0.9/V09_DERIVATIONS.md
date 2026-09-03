# v0.9 Detailed Derivations

This document expands the algebra used by the convergence and basis-management layer.

---

## A. Equal-width Gaussian saddle point

Consider the real exponent of $g_i^*g_j$:

$$
R(q)
=
-\frac12(q-q_i)^TA(q-q_i)
-\frac12(q-q_j)^TA(q-q_j).
$$

Then

$$
\nabla R
=
-A(q-q_i)-A(q-q_j).
$$

At the stationary point,

$$
A(2q-q_i-q_j)=0.
$$

For positive-definite $A$,

$$
\boxed{q_c=(q_i+q_j)/2.}
$$

The Hessian is

$$
\nabla\nabla R=-2A,
$$

which is negative definite.  Thus $q_c$ is the unique maximum of the product
magnitude.

---

## B. Complex first moment

For equal widths,

$$
g_i^*g_j
\propto
\exp\left[-(q-q_c)^TA(q-q_c)+i\Delta p^Tq\right],
$$

where

$$
\Delta p=p_j-p_i.
$$

Complete the square:

$$
-(q-q_c)^TA(q-q_c)+i\Delta p^T(q-q_c)
$$

$$
=
-(q-\mu)^TA(q-\mu)
-\frac14\Delta p^TA^{-1}\Delta p,
$$

with

$$
\boxed{\mu=q_c+\frac{i}{2}A^{-1}\Delta p.}
$$

Therefore

$$
\frac{\langle g_i|q|g_j\rangle}{S_{ij}}=\mu.
$$

---

## C. SPA0

Expand

$$
F(q)=F_c+\mathcal O(q-q_c).
$$

Then

$$
\langle g_i|F|g_j\rangle
\approx
F_c\langle g_i|g_j\rangle.
$$

Hence

$$
\boxed{M_{ij}^{(0)}=F_cS_{ij}.}
$$

---

## D. SPA1

Retain first order:

$$
F(q)
\approx
F_c
+
\sum_\alpha F_{c,\alpha}(q_\alpha-q_{c,\alpha}).
$$

Therefore

$$
M_{ij}^{(1)}
=
F_cS_{ij}
+
\sum_\alpha
F_{c,\alpha}
\langle g_i|(q_\alpha-q_{c,\alpha})|g_j\rangle.
$$

Use the first moment:

$$
\langle g_i|(q-q_c)|g_j\rangle
=(\mu-q_c)S_{ij}.
$$

Thus

$$
\boxed{
M_{ij}^{(1)}
=
\left[F_c+\nabla F_c\cdot(\mu-q_c)\right]S_{ij}.
}
$$

For

$$
F(q)=a+b^Tq,
$$

the Taylor expansion terminates exactly, so the SPA1 expression equals the exact
Gaussian integral.

---

## E. Second-order remainder

Let

$$
F(q)=F_c+g_c^T\delta q+\frac12\delta q^TH_F(\xi)\delta q,
$$

for some point $\xi$ in the Taylor remainder.

Then

$$
R_{ij}^{(2)}
=
\frac12
\langle g_i|
\delta q^TH_F(\xi)\delta q
|g_j\rangle.
$$

If the Hessian norm is bounded by $L$ over the important overlap region,

$$
|R_{ij}^{(2)}|
\lesssim
\frac{L}{2}
\langle g_i|\|\delta q\|^2|g_j\rangle.
$$

This makes explicit why broader Gaussians and rapidly varying electronic quantities
reduce the validity of low-order saddle-point approximations.

---

## F. Hermiticity of SPA1

For a Hermitian operator field,

$$
F(q)=F^\dagger(q),
$$

and a symmetric pair saddle,

$$
F_{ji}^{(0)}=(F_{ij}^{(0)})^*,
$$

$$
F_{ji,\alpha}^{(1)}=(F_{ij,\alpha}^{(1)})^*.
$$

Also

$$
S_{ji}=S_{ij}^*,
\qquad
\mu_{ji}-q_c=(\mu_{ij}-q_c)^*.
$$

Therefore

$$
M_{ji}^{(1)}=(M_{ij}^{(1)})^*.
$$

---

## G. Projection after pruning

Let retained basis vectors be columns of $\Phi_K$.  Minimize

$$
\mathcal L
=
\|\Phi C-\Phi_KC'\|^2.
$$

Differentiate with respect to $(C')^*$:

$$
\frac{\partial\mathcal L}{\partial(C')^*}
=
-S_{K,\mathrm{all}}C+S_{KK}C'=0.
$$

Hence

$$
\boxed{S_{KK}C'=S_{K,\mathrm{all}}C.}
$$

The projected norm is

$$
N_K=(C')^\dagger S_{KK}C'.
$$

Orthogonal projection cannot increase the norm, so

$$
\boxed{N_{\mathrm{old}}-N_K\ge0.}
$$

The equality holds for exact redundancy.

---

## H. Smallest overlap eigenvector and redundancy

Suppose

$$
Su=\lambda u
$$

with $\lambda\ll1$.

The norm of the linear combination

$$
|\Xi\rangle=\sum_i u_i|G_i\rangle
$$

is

$$
\langle\Xi|\Xi\rangle
=u^\dagger Su
=\lambda u^\dagger u.
$$

For normalized $u$,

$$
\boxed{\|\Xi\|^2=\lambda.}
$$

Thus a tiny overlap eigenvalue directly identifies an almost-null linear combination.

---

## I. Canonical orthogonalization

Let

$$
S=U\Lambda U^\dagger.
$$

Define

$$
X=U\Lambda^{-1/2}.
$$

Then

$$
X^\dagger SX
=
\Lambda^{-1/2}U^\dagger U\Lambda U^\dagger U\Lambda^{-1/2}
=I.
$$

After discarding small eigenvalues, the same identity holds in the retained subspace.

---

## J. Integrated coupling action

Consider two states and neglect all terms except one nonadiabatic source term:

$$
\dot c_b=-\eta_{ba}(t)c_a,
$$

where

$$
\eta_{ba}=\dot q\cdot d_{ba}.
$$

If $c_a$ varies slowly over a short coupling interval,

$$
\Delta c_b
\approx
-c_a\int\eta_{ba}(t)dt.
$$

Taking absolute values motivates

$$
\boxed{
\mathcal A_{ba}
=
\int|\eta_{ba}(t)|dt.
}
$$

This is not a transition probability.  It is a dimensionless measure of cumulative
first-order coupling exposure.

---

## K. Time-step invariance of the action sum

For constant coupling rate $\eta$,

$$
\mathcal A_N
=
\sum_{n=1}^N|\eta|\Delta t
=N|\eta|\Delta t.
$$

Since

$$
t=N\Delta t,
$$

$$
\boxed{\mathcal A(t)=|\eta|t.}
$$

Thus the ideal trigger time

$$
\mathcal A(t_*)=\mathcal A_{\mathrm{spawn}}
$$

is

$$
\boxed{t_*=\mathcal A_{\mathrm{spawn}}/|\eta|,}
$$

independent of $\Delta t$, up to discrete crossing of the threshold.

---

## L. Observed convergence order

Assume

$$
Q(h)=Q^*+Ch^p.
$$

Then

$$
Q(h)-Q(h/r)
=Ch^p(1-r^{-p}),
$$

and

$$
Q(h/r)-Q(h/r^2)
=Ch^pr^{-p}(1-r^{-p}).
$$

Taking the ratio,

$$
\frac{|Q(h)-Q(h/r)|}{|Q(h/r)-Q(h/r^2)|}=r^p.
$$

Therefore

$$
\boxed{
p=\frac{\ln(e_h/e_{h/r})}{\ln r}.}
$$

---

## M. Population reference error

For two-state population vectors

$$
P^{A}=(P_0^A,P_1^A),
\qquad
P^{Q}=(P_0^Q,P_1^Q),
$$

define

$$
\boxed{
\epsilon_P
=
\sqrt{
(P_0^A-P_0^Q)^2
+
(P_1^A-P_1^Q)^2
}.
}
$$

This is the v0.9 default scalar population error against the exact two-dimensional
reference.
