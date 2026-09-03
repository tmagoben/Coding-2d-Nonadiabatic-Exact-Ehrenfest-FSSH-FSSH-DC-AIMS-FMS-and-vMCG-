# v0.7 Detailed Derivations

## A. Polar link gauge covariance

Let

$$
O=U\Sigma V^\dagger.
$$

The unitary polar factor is

$$
L=UV^\dagger.
$$

Apply unitary node gauges $G_a,G_b$:

$$
O'=G_a^\dagger O G_b.
$$

Because left and right multiplication by unitary matrices does not change the singular
values, one valid SVD is

$$
O'
=
(G_a^\dagger U)
\Sigma
(G_b^\dagger V)^\dagger.
$$

Therefore

$$
L'
=
G_a^\dagger U V^\dagger G_b,
$$

so

$$
\boxed{
L'=G_a^\dagger L G_b.
}
$$

---

## B. Wilson-loop transformation

For a three-node loop,

$$
W=U_{01}U_{12}U_{20}.
$$

After local gauge transformations,

$$
W'
=
(G_0^\dagger U_{01}G_1)
(G_1^\dagger U_{12}G_2)
(G_2^\dagger U_{20}G_0).
$$

Adjacent gauge matrices cancel:

$$
G_1G_1^\dagger=I,
\qquad
G_2G_2^\dagger=I.
$$

Hence

$$
\boxed{
W'=G_0^\dagger W G_0.
}
$$

Similarity transformations preserve eigenvalues and trace.

---

## C. Spanning-tree gauge

We require a transformed tree link to be the identity:

$$
G_p^\dagger U_{pc}G_c=I.
$$

Multiply from the left by $G_p$:

$$
U_{pc}G_c=G_p.
$$

Then

$$
\boxed{
G_c=U_{pc}^\dagger G_p.
}
$$

Recursive application along the tree determines every node gauge from the root.

---

## D. Synchronization objective

For one edge,

$$
\|G_u^\dagger U_{uv}G_v-I\|_F^2.
$$

Because the matrices are unitary,

$$
\|A-I\|_F^2
=
2m-2\operatorname{ReTr}A.
$$

Therefore minimizing the total objective is equivalent to maximizing

$$
\sum_{(u,v)}
w_{uv}
\operatorname{ReTr}
(G_u^\dagger U_{uv}G_v).
$$

Hold every gauge except $G_u$ fixed.  Define

$$
A_u
=
\sum_v
w_{uv}U_{uv}G_v.
$$

The node-dependent objective becomes

$$
\max_{G_u\in U(m)}
\operatorname{ReTr}(G_u^\dagger A_u).
$$

If

$$
A_u=X\Sigma Y^\dagger,
$$

the Procrustes solution is

$$
\boxed{
G_u=XY^\dagger.
}
$$

That is precisely the nearest-unitary/polar factor used in the code.

---

## E. Coefficient transport covariance

Suppose

$$
c_u=U_{uv}c_v.
$$

Under local gauges,

$$
c_v'=G_v^\dagger c_v
$$

and

$$
U_{uv}'=G_u^\dagger U_{uv}G_v.
$$

Then

$$
c_u'
=U_{uv}'c_v'
$$

$$
=
G_u^\dagger U_{uv}G_vG_v^\dagger c_v
$$

$$
=
G_u^\dagger U_{uv}c_v.
$$

Therefore

$$
\boxed{
c_u'=G_u^\dagger c_u.
}
$$

The transported physical state is gauge covariant.

---

## F. Derivative-Hamiltonian matrix

Differentiate

$$
H|\phi_j\rangle=E_j|\phi_j\rangle.
$$

Project with $\langle\phi_i|$:

$$
\langle\phi_i|\partial H|\phi_j\rangle
+
E_i\langle\phi_i|\partial\phi_j\rangle
=
\delta_{ij}\partial E_j
+
E_j\langle\phi_i|\partial\phi_j\rangle.
$$

For $i=j$,

$$
\boxed{
F_{ii}=\partial E_i.
}
$$

For $i\ne j$,

$$
\boxed{
F_{ij}
=(E_j-E_i)d_{ij}.
}
$$

Thus $F$ remains well defined as an electronic operator matrix even when one later
rotates the electronic frame.

---

## G. Gauge transformation of the derivative connection

Let

$$
\Phi'=\Phi G(q).
$$

Then

$$
d'
=
\Phi'^\dagger\partial\Phi'.
$$

Substitute:

$$
d'
=
G^\dagger\Phi^\dagger
\partial(\Phi G).
$$

Use the product rule:

$$
\partial(\Phi G)
=(\partial\Phi)G+\Phi(\partial G).
$$

Hence

$$
d'
=
G^\dagger dG
+
G^\dagger(\partial G),
$$

so

$$
\boxed{
d'=G^\dagger dG+G^\dagger\partial G.
}
$$

This is why v0.7 treats overlap links as the discrete connection rather than rotating
NACs as ordinary operators.

---

## H. Gauge invariance of pair electronic overlap

Transport two local coefficient vectors to reference node $r$:

$$
\tilde c_i=T_{i\to r}c_i,
\qquad
\tilde c_j=T_{j\to r}c_j.
$$

Under local gauge transformations,

$$
\tilde c_i' = G_r^\dagger\tilde c_i,
\qquad
\tilde c_j' = G_r^\dagger\tilde c_j.
$$

Therefore

$$
(\tilde c_i')^\dagger\tilde c_j'
=
\tilde c_i^\dagger
G_rG_r^\dagger
\tilde c_j,
$$

so

$$
\boxed{
s_{ij}^{(e)'}=s_{ij}^{(e)}.
}
$$

---

## I. Gauge invariance of pair potential factor

At the reference node,

$$
H_r'=G_r^\dagger H_rG_r.
$$

Thus

$$
(\tilde c_i')^\dagger
H_r'
\tilde c_j'
$$

$$
=
\tilde c_i^\dagger
G_r
G_r^\dagger H_rG_r
G_r^\dagger
\tilde c_j
$$

$$
=
\boxed{
\tilde c_i^\dagger H_r\tilde c_j.
}
$$

The same proof applies to every derivative-Hamiltonian matrix $F_\alpha$.

---

## J. Gauge invariance of the graph-Gaussian pair matrix

The nuclear Gaussian quantities

$$
S_{ij}^{(n)},
\qquad
T_{ij}^{(n)}
$$

are electronic-gauge independent.

The electronic factors

$$
s_{ij}^{(e)},
\qquad
v_{ij}^{(e)}
$$

were proved invariant above.

Therefore

$$
S_{ij}
=S_{ij}^{(n)}s_{ij}^{(e)}
$$

and

$$
H_{ij}
=T_{ij}^{(n)}s_{ij}^{(e)}
+S_{ij}^{(n)}v_{ij}^{(e)}
$$

are invariant under arbitrary local node gauges.

This is tested with independent random $U(2)$ transformations at every graph node.

---

## K. Cayley norm conservation

The generalized Crank-Nicolson equation is

$$
\left(S+\frac{i\Delta t}{2}H\right)C_{n+1}
=
\left(S-\frac{i\Delta t}{2}H\right)C_n.
$$

For Hermitian $S$ and $H$, the corresponding propagator is $S$-unitary:

$$
U^\dagger S U=S.
$$

Hence

$$
\boxed{
C_{n+1}^\dagger S C_{n+1}
=
C_n^\dagger S C_n.
}
$$

The v0.7 static branched-basis regression verifies this to roundoff.
