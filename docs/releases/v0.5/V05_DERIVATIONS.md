# v0.5 Detailed Derivations

---

## A. Generalized-coordinate mass matrix

Let

$$
R=R_0+Jq.
$$

Then

$$
\dot R=J\dot q.
$$

The Cartesian kinetic energy is

$$
T
=
\frac12
\dot R^TM_R\dot R.
$$

Substitute:

$$
T
=
\frac12
\dot q^TJ^TM_RJ\dot q.
$$

Hence

$$
\boxed{
M_q=J^TM_RJ.
}
$$

The generalized canonical momentum is

$$
p=\frac{\partial T}{\partial\dot q}=M_q\dot q.
$$

Therefore

$$
\boxed{
\dot q=M_q^{-1}p.
}
$$

---

## B. Projection of gradients

For a scalar electronic energy,

$$
E=E(R(q)).
$$

By the chain rule,

$$
\frac{\partial E}{\partial q_\alpha}
=
\sum_A
\frac{\partial E}{\partial R_A}
\frac{\partial R_A}{\partial q_\alpha}.
$$

In matrix notation,

$$
\boxed{
g_q=J^Tg_R.
}
$$

The same coordinate transformation applies to the derivative-coupling covector:

$$
\boxed{
d_q=J^Td_R.
}
$$

---

## C. Complex overlap centroid

For equal-width Gaussians,

$$
g_i^*g_j
\propto
\exp
\left[
-\frac12(q-q_i)^TA(q-q_i)
-\frac12(q-q_j)^TA(q-q_j)
+
i(p_j-p_i)^Tq
\right].
$$

Let

$$
\bar q=\frac{q_i+q_j}{2},
\qquad
\Delta p=p_j-p_i.
$$

The quadratic position part is

$$
-\frac12
[(q-q_i)^TA(q-q_i)+(q-q_j)^TA(q-q_j)]
$$

$$
=
-(q-\bar q)^TA(q-\bar q)
-\frac14(q_i-q_j)^TA(q_i-q_j).
$$

The remaining linear phase is

$$
i\Delta p^Tq.
$$

Complete the square:

$$
-(q-\bar q)^TA(q-\bar q)
+
i\Delta p^T(q-\bar q)
$$

$$
=
-(q-\mu)^TA(q-\mu)
-
\frac14\Delta p^TA^{-1}\Delta p,
$$

where

$$
\boxed{
\mu
=
\bar q
+
\frac{i}{2}A^{-1}\Delta p.
}
$$

This is a complex centroid because the cross density $g_i^*g_j$ is not a probability
density.

---

## D. Gradient matrix element

For Gaussian $j$,

$$
\nabla g_j
=
[-A(q-q_j)+ip_j]g_j.
$$

Because the normalized cross-Gaussian first moment is $\mu$,

$$
\langle q\rangle_{ij}=\mu,
$$

so

$$
\boxed{
\langle g_i|\nabla g_j\rangle
=
[-A(\mu-q_j)+ip_j]S_{ij}.
}
$$

---

## E. Kinetic matrix element

Write

$$
B=M_q^{-1}.
$$

By integration by parts,

$$
T_{ij}
=
\frac12
\int
(\nabla g_i)^TB(\nabla g_j)
\,dq.
$$

Define

$$
u_i(q)
=
-A(q-q_i)-ip_i,
$$

$$
u_j(q)
=
-A(q-q_j)+ip_j.
$$

Then

$$
T_{ij}
=
\frac12S_{ij}
\langle
u_i^TBu_j
\rangle_{ij}.
$$

Decompose

$$
q=\mu+\delta.
$$

The cross Gaussian has covariance

$$
\boxed{
\langle\delta\delta^T\rangle_{ij}
=
\frac12A^{-1}.
}
$$

Therefore

$$
\langle
u_i^TBu_j
\rangle
=
u_i(\mu)^TBu_j(\mu)
+
\operatorname{Tr}
\left[
A B A
\frac12A^{-1}
\right].
$$

Since

$$
ABA A^{-1}=AB,
$$

the trace contribution is

$$
\frac12\operatorname{Tr}(BA).
$$

Hence

$$
\boxed{
T_{ij}
=
\frac12S_{ij}
\left[
u_i(\mu)^TBu_j(\mu)
+
\frac12\operatorname{Tr}(BA)
\right].
}
$$

---

## F. First-order NAC Hamiltonian term

The locally constant first-order adiabatic kinetic coupling is

$$
\hat H^{(1)}_{ab}
=
-
d_{ab}^TB\nabla.
$$

Therefore

$$
H^{(1)}_{ia,jb}
=
-
d_{ab}^TB
\langle g_i|\nabla g_j\rangle.
$$

Using the gradient matrix element,

$$
\boxed{
H^{(1)}_{ia,jb}
=
-
d_{ab}^TB
G_{ij}.
}
$$

Now use

$$
d_{ba}=-d_{ab}
$$

and

$$
G_{ji}^*=-G_{ij}.
$$

Then

$$
(H^{(1)}_{jb,ia})^*
=
-
(-d_{ab})^TB(-G_{ij})
=
H^{(1)}_{ia,jb}.
$$

Thus the first-order local coupling is Hermitian.

---

## G. Locally constant $d^2$ contribution

For nuclear dimensions $\alpha,\beta$,

$$
D^{(2)}_{ab}
=
\sum_{c,\alpha,\beta}
d_{ac,\alpha}
B_{\alpha\beta}
d_{cb,\beta}.
$$

In matrix notation, if $d_\alpha$ is the electronic coupling matrix associated with
coordinate $\alpha$,

$$
\boxed{
D^{(2)}
=
\sum_{\alpha\beta}
B_{\alpha\beta}
d_\alpha d_\beta.
}
$$

The contribution to the Hamiltonian is

$$
\boxed{
H^{(2)}_{ia,jb}
=
-\frac12
D^{(2)}_{ab}
S_{ij}.
}
$$

For a two-state real problem and diagonal scalar mass, the diagonal term becomes the
usual positive diagonal Born-Oppenheimer-like contribution because $d^2$ is negative.

---

## H. Moving-basis matrix

For fixed $A$,

$$
g_j
=
N
e^{-\frac12(q-q_j)^TA(q-q_j)+ip_j^T(q-q_j)}.
$$

Differentiate with respect to time:

$$
\dot g_j
=
\left[
(A(q-q_j)-ip_j)^T\dot q_j
+
i(q-q_j)^T\dot p_j
\right]g_j.
$$

Take the matrix element with $g_i$ and replace $q$ by the cross centroid $\mu$:

$$
\boxed{
T^{\mathrm{basis}}_{ij}
=
S_{ij}
\left[
(A(\mu-q_j)-ip_j)^T\dot q_j
+
i(\mu-q_j)^T\dot p_j
\right].
}
$$

---

## I. General-mass NAC momentum rescaling

Choose a momentum adjustment

$$
p'=p+\lambda n.
$$

Let

$$
B=M^{-1}.
$$

Energy conservation gives

$$
\frac12(p+\lambda n)^TB(p+\lambda n)+E_b
=
\frac12p^TBp+E_a.
$$

Expand:

$$
\frac12
\left[
2\lambda p^TBn
+
\lambda^2n^TBn
\right]
+
E_b-E_a
=
0.
$$

Multiply by two:

$$
(n^TBn)\lambda^2
+
2(p^TBn)\lambda
+
2(E_b-E_a)
=
0.
$$

Therefore,

$$
\boxed{
\lambda
=
\frac{
-p^TBn
\pm
\sqrt{
(p^TBn)^2
-
2(n^TBn)(E_b-E_a)
}
}{
n^TBn
}.
}
$$

---

## J. Why a local electronic approximation is necessary

A formally exact Gaussian pair matrix element requires integrals such as

$$
\int
g_i^*(q)
E_a(q)
g_j(q)
dq
$$

and

$$
\int
g_i^*(q)
d_{ab}(q)\cdot\nabla g_j(q)
dq.
$$

For an on-the-fly ab initio method, evaluating $E$ and $d$ on all points in the
Gaussian overlap region is computationally prohibitive.

The centroid approximation replaces

$$
E(q)\rightarrow E(q_{ij}),
$$

$$
d(q)\rightarrow d(q_{ij}),
$$

so the remaining nuclear integrals are analytic.

The approximation is controlled by how slowly the electronic quantities vary over the
Gaussian overlap region. It should be checked by narrowing widths, comparing selected
pair matrix elements against additional backend points, or moving to a higher-order
saddle-point approximation in a later release.
