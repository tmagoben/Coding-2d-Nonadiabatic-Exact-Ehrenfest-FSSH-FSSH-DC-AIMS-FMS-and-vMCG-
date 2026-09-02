# v0.7 Theory: Graph-Based Electronic Gauge Transport for Branched Gaussian Dynamics

Version 0.7 addresses a limitation that appears as soon as more than one Gaussian
trajectory basis function (TBF) is propagated at once.

A single trajectory has a natural sequence

```text
R0 -> R1 -> R2 -> R3
```

and v0.6 can track its electronic states sequentially.  A spawned Gaussian basis is
not a line.  It is a network:

```text
        TBF 1
       /     \
 centroid   centroid
    /           \
 TBF 0 -------- TBF 2
```

Different TBF centers and Gaussian-pair centroids may all require electronic
wavefunctions.  There is then no unique "previous geometry" that can define the
phase/gauge everywhere.

v0.7 replaces the single tracking chain by a **discrete electronic gauge graph**.

Atomic units are used for the dynamics.

---

# 1. Electronic frame at each node

At graph node $n$, let

$$
\boxed{
\Phi_n
=
\left(
|\phi_{n1}\rangle,
\ldots,
|\phi_{nm}\rangle
\right)
}
$$

contain an orthonormal set of $m$ electronic states as columns.

The frame is not unique.  A local unitary change of basis gives

$$
\boxed{
\Phi_n\rightarrow\Phi_n G_n,
\qquad
G_n\in U(m).
}
$$

For isolated nondegenerate real states, this reduces to independent signs.
For a degenerate or near-degenerate manifold, the full unitary freedom matters.

---

# 2. Edge overlap matrix

For neighboring nodes $u$ and $v$, define

$$
\boxed{
O_{uv}
=
\Phi_u^\dagger\Phi_v.
}
$$

In an ab initio calculation, this is a many-electron overlap matrix.  v0.6 already
constructs it from PySCF SA-CASSCF wavefunctions in nonorthogonal orbital bases.

Under local gauge transformations,

$$
O_{uv}
\rightarrow
G_u^\dagger O_{uv}G_v.
$$

If the selected state subspaces at the two nodes are identical, $O_{uv}$ is unitary.
For finite geometry steps or an incomplete state manifold it need not be exactly
unitary.

---

# 3. Polar/unitary link

Take the singular-value decomposition

$$
O_{uv}=A\Sigma B^\dagger.
$$

The nearest unitary matrix in Frobenius norm is

$$
\boxed{
U_{uv}=AB^\dagger.
}
$$

This is the unitary polar factor of the overlap matrix.

It inherits the same gauge transformation law:

$$
\boxed{
U_{uv}
\rightarrow
G_u^\dagger U_{uv}G_v.
}
$$

Therefore $U_{uv}$ is the discrete analogue of an electronic gauge connection.

The reverse link is

$$
\boxed{
U_{vu}=U_{uv}^\dagger.
}
$$

---

# 4. Transport of electronic coefficients

Suppose an electronic state is represented at node $v$ by coefficients $c_v$:

$$
|\Psi\rangle=\Phi_v c_v.
$$

Because

$$
\Phi_v\approx\Phi_uU_{uv},
$$

its coefficients expressed in node-$u$ coordinates are

$$
\boxed{
c_u=U_{uv}c_v.
}
$$

Along a path

$$
n_0\rightarrow n_1\rightarrow\cdots\rightarrow n_k,
$$

coefficients are transported by the ordered product

$$
\boxed{
T_{n_0\rightarrow n_k}
=
U_{n_k n_{k-1}}
\cdots
U_{n_2n_1}
U_{n_1n_0}.
}
$$

This is implemented directly in `gauge_graph.py`.

---

# 5. Wilson loop / holonomy

For a closed cycle

$$
\mathcal C:
0\rightarrow1\rightarrow\cdots\rightarrow N\rightarrow0,
$$

define the Wilson product

$$
\boxed{
W_{\mathcal C}
=
U_{01}U_{12}\cdots U_{N0}.
}
$$

Under node-local gauges,

$$
W_{\mathcal C}
\rightarrow
G_0^\dagger W_{\mathcal C}G_0.
$$

Therefore the following are gauge invariant:

$$
\boxed{
\operatorname{Tr}W_{\mathcal C}
}
$$

and the eigenvalues of $W_{\mathcal C}$.

For a single isolated state, $m=1$, so

$$
W_{\mathcal C}=e^{i\gamma_{\mathcal C}},
$$

which is the discrete Berry phase.

For a degenerate manifold, $W_{\mathcal C}$ is matrix valued and represents a
non-Abelian holonomy.

---

# 6. Conical-intersection loop

For the real lower adiabatic state around the two-dimensional CI introduced in v0.4,
the Berry phase is

$$
\gamma=\pi.
$$

Therefore

$$
\boxed{
W_{\mathcal C}=-1.
}
$$

v0.7 constructs a ring graph from the analytic CI eigenvectors and verifies this
directly.

This result is important conceptually:

> a globally consistent gauge does **not** mean forcing every loop product to the
> identity.

A nontrivial physical holonomy must remain visible.

---

# 7. Spanning-tree gauge

Choose a graph root $r$ and a spanning tree.

Set

$$
G_r=I.
$$

For a tree edge from parent $p$ to child $c$, choose $G_c$ so that the transformed
link is the identity:

$$
G_p^\dagger U_{pc}G_c=I.
$$

Solving gives

$$
\boxed{
G_c=U_{pc}^\dagger G_p.
}
$$

Thus every tree edge can be made exactly trivial.

Any remaining non-tree edge closes a cycle.  Its residual transformed link is then a
representation of the cycle holonomy.

This cleanly separates

```text
removable local gauge variation
```

from

```text
non-removable loop holonomy.
```

---

# 8. Why a tree alone is not enough

A spanning tree contains no closed loops.

Consequently one can always make all of its links equal to $I$, even if the full
graph contains a Berry phase or non-Abelian curvature.

The physically interesting information therefore lives on the **chords** that close
fundamental cycles.

v0.7 constructs a fundamental cycle basis relative to a selected spanning tree and
computes the Wilson matrix for each cycle.

---

# 9. Global gauge synchronization

For noisy finite-step overlaps, it can be undesirable to concentrate all graph
inconsistency on one non-tree edge.

v0.7 also solves the anchored unitary synchronization problem

$$
\boxed{
\min_{\{G_n\}}
\sum_{(u,v)}
w_{uv}
\left\|
G_u^\dagger U_{uv}G_v-I
\right\|_F^2.
}
$$

The root gauge is fixed,

$$
G_r=I,
$$

to remove the arbitrary global unitary freedom.

This objective attempts to make neighboring frames as mutually smooth as possible.

If every Wilson loop is trivial, the minimum can be zero.

If the graph carries physical holonomy, the minimum remains nonzero.

The algorithm therefore **distributes** unavoidable loop frustration rather than
erasing it.

---

# 10. Coordinate-descent synchronization update

Holding all neighboring gauges fixed, the terms involving node $u$ reduce to
maximizing

$$
\operatorname{ReTr}(G_u^\dagger A_u),
$$

where

$$
\boxed{
A_u
=
\sum_{v\in\mathcal N(u)}
w_{uv}U_{uv}G_v.
}
$$

The maximizing unitary is the polar factor of $A_u$:

$$
\boxed{
G_u
=
\operatorname{polar}(A_u).
}
$$

v0.7 uses this block update as the refinement step.  To avoid a poor stationary point
from one arbitrary tree gauge, the implementation compares a spectral connection-matrix
initialization, the spanning-tree initialization, and a small deterministic set of seeded
random-unitary starts.  Each candidate is refined by the same transparent block update and
the lowest-objective result is retained.

The purpose is pedagogical transparency and robust small-graph behavior, not large-scale
graph optimization.

---

# 11. Operator matrices rather than raw state labels

At a node, the adiabatic electronic Hamiltonian is

$$
\boxed{
H_e=\operatorname{diag}(E_1,\ldots,E_m).
}
$$

A unitary local frame transformation gives

$$
\boxed{
H_e' = G^\dagger H_eG.
}
$$

After a general rotation inside a near-degenerate manifold, the Hamiltonian need no
longer be diagonal.

This is expected.  The local frame has become diabatic/quasi-diabatic rather than
adiabatic.

---

# 12. Derivative-Hamiltonian operator

In an adiabatic basis,

$$
F_{ij}^{(\alpha)}
=
\left\langle
\phi_i
\middle|
\frac{\partial H_e}{\partial q_\alpha}
\middle|
\phi_j
\right\rangle.
$$

For the diagonal elements,

$$
\boxed{
F_{ii}^{(\alpha)}
=
\frac{\partial E_i}{\partial q_\alpha}.
}
$$

For $i\ne j$,

$$
\boxed{
F_{ij}^{(\alpha)}
=
(E_j-E_i)d_{ij}^{(\alpha)}.
}
$$

Thus v0.7 reconstructs the full operator matrix from the energies, adiabatic gradients,
and NACs already returned by the v0.5/v0.6 provider layer.

Under a node-local unitary transformation,

$$
\boxed{
F_\alpha' = G^\dagger F_\alpha G.
}
$$

Unlike the derivative-coupling connection itself, this is an ordinary covariant
operator transformation.

---

# 13. Why the NAC cannot simply be rotated as an operator

The derivative coupling is a gauge connection:

$$
\boxed{
d_\alpha
=
\Phi^\dagger\partial_\alpha\Phi.
}
$$

For a coordinate-dependent gauge $G(q)$,

$$
\boxed{
d_\alpha'
=
G^\dagger d_\alpha G
+
G^\dagger\partial_\alpha G.
}
$$

The inhomogeneous second term is essential.

Therefore v0.7 does **not** claim that an arbitrarily graph-rotated NAC is obtained by
$G^\dagger dG$ alone.

Instead:

- overlap links represent the discrete connection;
- $H_e$ and $F_\alpha$ are rotated as ordinary electronic operators;
- Wilson loops diagnose connection holonomy.

This distinction is central to the scientific correctness of the release.

---

# 14. TBF-center and pair-centroid graph

For Gaussian TBFs $i$ and $j$, a local pair approximation may evaluate electronic
structure at their midpoint/centroid node $c_{ij}$.

A natural graph therefore contains

```text
TBF-center nodes
+
pair-centroid nodes
```

with edges

$$
i\leftrightarrow c_{ij},
\qquad
j\leftrightarrow c_{ij}.
$$

For three TBFs, the pair-centroid graph automatically contains loops.

Those loops make gauge consistency a network problem rather than a sequential-path
problem.

`pyscf_gauge_graph.py` includes a helper that constructs exactly this connectivity.

---

# 15. Transport to one pair reference frame

Suppose TBF $i$ carries electronic coefficients $c_i$ in its local node frame and TBF
$j$ carries $c_j$.

Transport both to centroid/reference node $c$:

$$
\boxed{
\tilde c_i
=
T_{i\rightarrow c}c_i,
\qquad
\tilde c_j
=
T_{j\rightarrow c}c_j.
}
$$

Then define the electronic overlap factor

$$
\boxed{
s_{ij}^{(e)}
=
\tilde c_i^\dagger\tilde c_j.
}
$$

The local electronic potential factor is

$$
\boxed{
v_{ij}^{(e)}
=
\tilde c_i^\dagger
H_e(c)
\tilde c_j.
}
$$

And the derivative-Hamiltonian factors are

$$
\boxed{
f_{ij,\alpha}^{(e)}
=
\tilde c_i^\dagger
F_\alpha(c)
\tilde c_j.
}
$$

These are invariant under arbitrary local unitary gauge changes when coefficients,
links, and operator matrices are transformed consistently.

That invariance is unit tested.

---

# 16. Discrete local-diabatic Gaussian pair approximation

Let the nuclear Gaussian overlap be

$$
S_{ij}^{(n)}
=
\langle g_i|g_j\rangle,
$$

and the nuclear kinetic matrix element be

$$
T_{ij}^{(n)}.
$$

v0.7 introduces the transparent approximation

$$
\boxed{
S_{ij}
\approx
S_{ij}^{(n)}
s_{ij}^{(e)}
}
$$

and

$$
\boxed{
H_{ij}
\approx
T_{ij}^{(n)}s_{ij}^{(e)}
+
S_{ij}^{(n)}v_{ij}^{(e)}.
}
$$

The electronic factors are evaluated only after both TBF electronic vectors have been
transported to one common graph node.

This makes the matrix element explicitly gauge covariant.

It is a **discrete local-diabatic approximation**, not the exact continuous AIMS
kinetic-coupling integral.

---

# 17. Why this is useful

The v0.5 centroid Hamiltonian used an adiabatic NAC expansion locally.

The v0.7 graph formulation supplies a complementary viewpoint:

```text
continuous derivative-coupling description
             versus
discrete wavefunction-overlap transport
```

The latter is particularly useful when analytical NACs are sharply peaked or when a
selected near-degenerate subspace is better handled as a locally diabatic block.

Local-diabatization ideas are well established in nonadiabatic trajectory methods;
v0.7 extends the conceptual structure from one trajectory to a graph of Gaussian
basis calculations.

---

# 18. Static graph-Gaussian propagation test

For a static nonorthogonal graph-Gaussian basis,

$$
iS\dot C=HC.
$$

v0.7 propagates one Cayley/Crank-Nicolson step as

$$
\boxed{
\left(S+\frac{i\Delta t}{2}H\right)C_{n+1}
=
\left(S-\frac{i\Delta t}{2}H\right)C_n.
}
$$

For Hermitian $H$ and positive $S$, this preserves

$$
\boxed{
C^\dagger SC
}
$$

up to linear-solver roundoff.

This provides a compact gauge-covariant branched-basis regression without pretending
that all production AIMS time-dependent matrix-element terms have already been solved.

---

# 19. Connection to PySCF

The graph can be built directly from v0.6
`CASSCFWavefunctionSnapshot` objects.

For each requested edge $(u,v)$,

$$
O_{uv}
=
\langle\Psi(u)|\Psi(v)\rangle
$$

is evaluated with the same nonorthogonal many-electron PySCF overlap machinery already
validated in v0.6.

Thus the intended direct-dynamics chain is

```text
PySCF SA-CASSCF snapshots
        |
many-electron overlaps
        |
unitary graph links
        |
cycle / holonomy diagnostics
        |
common-reference electronic transport
        |
gauge-covariant Gaussian pair matrices
```

---

# 20. What v0.7 does not claim

v0.7 does not claim that one can globally diabaticize an arbitrary molecular problem.

Nonzero curvature and CI topology can obstruct a globally flat gauge.

It also does not yet implement:

1. moving graph topology during a large production trajectory ensemble;
2. adaptive graph-edge selection from overlap/error estimates;
3. full time-dependent graph-Gaussian $T_{ij}=\langle G_i|\dot G_j\rangle$ including
   graph-link time derivatives;
4. higher-order saddle-point AIMS matrix elements;
5. a distributed graph optimizer for thousands of electronic nodes;
6. a non-Abelian graph connection generated directly from fully degenerate ab initio
   manifolds over a large molecular configuration space.

The v0.7 objective is more precise:

$$
\boxed{
\text{make gauge consistency a graph property and preserve physical holonomy.}
}
$$
