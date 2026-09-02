# v0.11 Detailed Derivations

This file expands the algebra used in `V11_THEORY.md`.

---

## A. Unequal-width Gaussian overlap

Let

$$
g_i(q)
=
N_i
e^{
-\frac12(q-q_i)^TA_i(q-q_i)
+
ip_i^T(q-q_i)
},
$$

and

$$
g_j(q)
=
N_j
e^{
-\frac12(q-q_j)^TA_j(q-q_j)
+
ip_j^T(q-q_j)
}.
$$

The complex-conjugated first Gaussian is

$$
g_i^*(q)
=
N_i
e^{
-\frac12(q-q_i)^TA_i(q-q_i)
-
ip_i^T(q-q_i)
}.
$$

Multiply:

$$
g_i^*g_j
=
N_iN_j
e^{
-\frac12q^TBq
+
\ell^Tq
+
c
},
$$

where

$$
\boxed{
B=A_i+A_j,
}
$$

$$
\boxed{
\ell=A_iq_i+A_jq_j+i(p_j-p_i),
}
$$

and

$$
\boxed{
c=
-\frac12q_i^TA_iq_i
-\frac12q_j^TA_jq_j
+
ip_i^Tq_i
-
ip_j^Tq_j.
}
$$

Complete the square:

$$
-\frac12q^TBq+\ell^Tq
=
-\frac12(q-\mu)^TB(q-\mu)
+
\frac12\ell^TB^{-1}\ell,
$$

with

$$
\boxed{
\mu=B^{-1}\ell.
}
$$

Use

$$
\int
e^{-\frac12x^TBx}
d^Dx
=
\frac{(2\pi)^{D/2}}{\sqrt{\det B}}.
$$

Therefore

$$
S_{ij}
=
N_iN_j
\frac{(2\pi)^{D/2}}{\sqrt{\det B}}
e^{
c+\frac12\ell^TB^{-1}\ell
}.
$$

Since

$$
N_iN_j
=
\frac{
(\det A_i\det A_j)^{1/4}
}{
\pi^{D/2}
},
$$

$$
\boxed{
S_{ij}
=
\frac{
2^{D/2}
(\det A_i\det A_j)^{1/4}
}{
\sqrt{\det(A_i+A_j)}
}
e^{
c+\frac12\ell^T(A_i+A_j)^{-1}\ell
}.
}
$$

---

## B. Equal-width reduction

Set

$$
A_i=A_j=A.
$$

Then

$$
B=2A.
$$

The prefactor becomes

$$
\frac{
2^{D/2}(\det A^2)^{1/4}
}{
\sqrt{\det(2A)}
}
=
\frac{
2^{D/2}(\det A)^{1/2}
}{
2^{D/2}(\det A)^{1/2}
}
=1.
$$

Also

$$
\mu
=
(2A)^{-1}
[A(q_i+q_j)+i(p_j-p_i)]
$$

so

$$
\boxed{
\mu
=
\frac{q_i+q_j}{2}
+
\frac{i}{2}A^{-1}(p_j-p_i).
}
$$

Thus the general formula exactly recovers the previous implementation.

---

## C. Cross covariance

After completing the square,

$$
g_i^*g_j
\propto
e^{-\frac12(q-\mu)^TB(q-\mu)}.
$$

The covariance of this Gaussian integral is

$$
\boxed{
\Sigma=B^{-1}.
}
$$

Thus

$$
\langle
(q-\mu)(q-\mu)^T
\rangle_{ij}
=
(A_i+A_j)^{-1}.
$$

---

## D. Real saddle point

The modulus ignores momentum phase:

$$
|g_i^*g_j|
\propto
e^{
-\frac12(q-q_i)^TA_i(q-q_i)
-\frac12(q-q_j)^TA_j(q-q_j)
}.
$$

Differentiate the exponent:

$$
A_i(q-q_i)+A_j(q-q_j)=0.
$$

Hence

$$
(A_i+A_j)q
=
A_iq_i+A_jq_j.
$$

Therefore

$$
\boxed{
q_s
=
(A_i+A_j)^{-1}
(A_iq_i+A_jq_j).
}
$$

---

## E. Gradient matrix element

Differentiate $g_j$:

$$
\nabla g_j
=
[-A_j(q-q_j)+ip_j]g_j.
$$

Then

$$
\langle g_i|\nabla g_j\rangle
=
S_{ij}
\left[
-A_j
(
\langle q\rangle_{ij}-q_j
)
+
ip_j
\right].
$$

Because

$$
\langle q\rangle_{ij}=\mu,
$$

$$
\boxed{
G_{ij}
=
[-A_j(\mu-q_j)+ip_j]S_{ij}.
}
$$

---

## F. Kinetic matrix element

For real symmetric mass metric,

$$
\hat T
=
-\frac12\nabla^TB_M\nabla,
\qquad
B_M=M^{-1}.
$$

Using vanishing boundary terms,

$$
T_{ij}
=
\frac12
\int
(\nabla g_i^*)^T
B_M
(\nabla g_j)dq.
$$

Define

$$
f_i(q)
=
-A_i(q-q_i)-ip_i,
$$

$$
f_j(q)
=
-A_j(q-q_j)+ip_j.
$$

Then

$$
T_{ij}
=
\frac12S_{ij}
\langle
f_i^TB_Mf_j
\rangle_{ij}.
$$

Write

$$
q=\mu+\delta.
$$

Then

$$
f_i=u_i-A_i\delta,
$$

$$
f_j=u_j-A_j\delta,
$$

where

$$
u_i=-A_i(\mu-q_i)-ip_i,
$$

$$
u_j=-A_j(\mu-q_j)+ip_j.
$$

The linear fluctuation terms vanish:

$$
\langle\delta\rangle=0.
$$

The quadratic term is

$$
\langle
\delta^T
A_iB_MA_j
\delta
\rangle
=
\operatorname{Tr}
[
A_iB_MA_j
\langle\delta\delta^T\rangle
].
$$

Using

$$
\langle\delta\delta^T\rangle=\Sigma,
$$

$$
\boxed{
T_{ij}
=
\frac12S_{ij}
[
u_i^TB_Mu_j
+
\operatorname{Tr}(A_iB_MA_j\Sigma)
].
}
$$

---

## G. Time derivative including width motion

The normalization is

$$
N_j
=
\left(
\frac{\det A_j}{\pi^D}
\right)^{1/4}.
$$

Therefore

$$
\frac{\dot N_j}{N_j}
=
\frac14
\frac{d}{dt}
\ln\det A_j.
$$

Use

$$
\frac{d}{dt}\ln\det A
=
\operatorname{Tr}(A^{-1}\dot A).
$$

Hence

$$
\boxed{
\frac{\dot N_j}{N_j}
=
\frac14
\operatorname{Tr}(A_j^{-1}\dot A_j).
}
$$

Now

$$
\xi_j=q-q_j.
$$

The width exponent is

$$
-\frac12\xi_j^TA_j\xi_j.
$$

Differentiate:

$$
\frac{d}{dt}
\left[
-\frac12\xi_j^TA_j\xi_j
\right]
=
\dot q_j^TA_j\xi_j
-
\frac12\xi_j^T\dot A_j\xi_j.
$$

The phase derivative is

$$
\frac{d}{dt}
[ip_j^T\xi_j]
=
i\dot p_j^T\xi_j
-
ip_j^T\dot q_j.
$$

Combine:

$$
\boxed{
\frac{\dot g_j}{g_j}
=
\frac14\operatorname{Tr}(A_j^{-1}\dot A_j)
+
(A_j\xi_j-ip_j)^T\dot q_j
+
i\xi_j^T\dot p_j
-
\frac12\xi_j^T\dot A_j\xi_j.
}
$$

For the cross integral,

$$
\langle\xi_j\rangle=\mu-q_j\equiv y.
$$

Also

$$
\langle
\xi_j^T\dot A_j\xi_j
\rangle
=
y^T\dot A_jy
+
\operatorname{Tr}(\dot A_j\Sigma).
$$

Substitute to obtain the implemented basis-time matrix element.

---

## H. Energy constraint for shifted spawn

Let

$$
p_c=p_p+\lambda n.
$$

The parent total classical energy is

$$
\mathcal E_p
=
E_a(q_p)
+
\frac12p_p^TBp_p.
$$

The child constraint is

$$
E_b(q_c)
+
\frac12
(p_p+\lambda n)^T
B
(p_p+\lambda n)
=
\mathcal E_p.
$$

Expand:

$$
E_b(q_c)
+
\frac12p_p^TBp_p
+
\lambda p_p^TBn
+
\frac12\lambda^2n^TBn
-
\mathcal E_p
=
0.
$$

Multiply by two:

$$
(n^TBn)\lambda^2
+
2(p_p^TBn)\lambda
+
p_p^TBp_p
+
2[E_b(q_c)-\mathcal E_p]
=
0.
$$

Thus

$$
a=n^TBn,
$$

$$
b=p_p^TBn,
$$

$$
c=p_p^TBp_p+2[E_b(q_c)-\mathcal E_p].
$$

The roots are

$$
\boxed{
\lambda
=
\frac{-b\pm\sqrt{b^2-ac}}{a}.
}
$$

---

## I. First-order coupling proxy

In the local adiabatic basis,

$$
H_e(q_s)
=
\operatorname{diag}(E_0,E_1,\ldots).
$$

For $a\ne b$,

$$
[H_e(q_s)]_{ab}=0.
$$

Expand:

$$
H_{ab}(q)
\approx
\sum_\alpha
F_{ab,\alpha}(q_s)
(q_\alpha-q_{s,\alpha}),
$$

where

$$
F_{ab,\alpha}
=
\langle\phi_a|
\partial_\alpha H_e
|\phi_b\rangle.
$$

Hellmann-Feynman differentiation gives

$$
\boxed{
F_{ab,\alpha}
=
(E_b-E_a)d_{ab,\alpha}.
}
$$

Therefore the Gaussian matrix element is approximately

$$
\langle g_p|
H_{ab}
|g_c\rangle
\approx
S_{pc}
\mathbf F_{ab}\cdot(\mu_{pc}-q_s).
$$

Hence the ranking proxy

$$
\boxed{
\mathcal V_{pc}
=
\left|
S_{pc}
\mathbf F_{ab}\cdot(\mu_{pc}-q_s)
\right|.
}
$$

---

## J. Novelty factor

Let

$$
s=\max_k|\langle g_k|g_c\rangle|
$$

over existing target-state TBFs.

The implemented novelty factor is

$$
\boxed{
\nu=(1-s^2)^\beta.
}
$$

Properties:

- if $s=0$, $\nu=1$;
- if $s\rightarrow1$, $\nu\rightarrow0$;
- $\nu$ decreases monotonically with redundancy.

The final ranking score is

$$
\boxed{
\mathcal J=\mathcal V_{pc}\nu.
}
$$

This novelty factor is a v0.11 heuristic, not part of the original optimal-spawning
paper.

---

## K. Canonical coefficient participation ratio

Diagonalize

$$
S=U\Lambda U^\dagger.
$$

Define orthonormal basis

$$
\chi=\Phi U\Lambda^{-1/2}.
$$

Then

$$
\Psi=\Phi C.
$$

Require

$$
\Phi C
=
\Phi U\Lambda^{-1/2}d.
$$

Therefore

$$
C=U\Lambda^{-1/2}d.
$$

Multiply by

$$
\Lambda^{1/2}U^\dagger:
$$

$$
\boxed{
d=\Lambda^{1/2}U^\dagger C.
}
$$

The norm is

$$
d^\dagger d
=
C^\dagger U\Lambda U^\dagger C
=
C^\dagger SC.
$$

Normalize

$$
p_k=\frac{|d_k|^2}{d^\dagger d}.
$$

Then

$$
\boxed{
N_{\rm part}
=
\frac{1}{\sum_kp_k^2}.
}
$$

This gives the effective number of occupied canonical directions in the
nonorthogonal basis.
