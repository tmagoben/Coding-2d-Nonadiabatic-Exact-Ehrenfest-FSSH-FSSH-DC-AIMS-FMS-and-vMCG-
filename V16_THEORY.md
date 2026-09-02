# v0.16 Theory: Persistent Locality Graphs, Sparse Propagation, and Provider-Aware Cost

v0.16 moves the Gaussian dynamics layer from a dense global pair representation toward
an explicitly local graph.

The release keeps the v0.15 analytic LVC physics and TDSE-defect controller, but changes
which Gaussian pairs are represented in the projected matrices.

The progression is

```text
v0.15
shared dense pair cache
+ incremental dense matrix updates
+ cost-aware residual growth

        ↓

v0.16
persistent overlap-locality graph
+ safe KD-tree geometric pre-screen
+ sparse S/H/T matrices
+ sparse Cayley solves
+ dense a-posteriori sparsification audit
+ local-degree cost model
+ electronic-structure cache cost term
```

## 1. Conservative overlap locality bound

For two normalized real-width Gaussian TBFs,

$$
g_i(R),\qquad g_j(R),
$$

define positive-definite width matrices

$$
A_i>0,\qquad A_j>0.
$$

The exact overlap magnitude is no larger than the zero-momentum spatial overlap.

The displacement term contains

$$
H_{ij}
=
(A_i^{-1}+A_j^{-1})^{-1}.
$$

If

$$
a_i=\lambda_{\min}(A_i),
\qquad
a_j=\lambda_{\min}(A_j),
$$

then

$$
\lambda_{\min}(H_{ij})
\ge
\frac{1}{1/a_i+1/a_j}.
$$

Therefore

$$
\boxed{
|S_{ij}|
\le
\exp
\left[
-\frac12h_{ij}
\|q_i-q_j\|^2
\right]
}
$$

with

$$
h_{ij}
=
\frac{1}{1/a_i+1/a_j}.
$$

Momentum mismatch can only reduce the exact overlap magnitude, so omitting it from the
screen is conservative.

## 2. Safe global KD-tree radius

Let

$$
a_{\min}
=
\min_i a_i.
$$

Then

$$
h_{ij}\ge a_{\min}/2.
$$

For an overlap threshold $\tau$,

$$
|S_{ij}|
\le
\exp
\left[
-\frac{a_{\min}}{4}
\|q_i-q_j\|^2
\right].
$$

Hence any pair satisfying

$$
\|q_i-q_j\|
>
\boxed{
R_\tau
=
\sqrt{
-4\ln\tau/a_{\min}
}
}
$$

cannot exceed the threshold.

v0.16 uses this radius in `scipy.spatial.cKDTree.query_pairs`.

Pairs outside $R_\tau$ are rejected without a pair solve.

Pairs inside the global radius receive the tighter pair-specific bound before an exact
overlap is evaluated.

## 3. Edge hysteresis

A single overlap cutoff can cause graph flicker as moving TBFs cross the threshold.

v0.16 uses

$$
\tau_{\rm enter}
>
\tau_{\rm exit}.
$$

The release values are

$$
\boxed{
\tau_{\rm enter}=0.03,
\qquad
\tau_{\rm exit}=0.015.
}
$$

A new edge must exceed $0.03$.

An existing edge remains until it falls below $0.015$.

Graph identity is keyed by persistent TBF `uid`, not transient list index.

## 4. Sparse projected matrices

Every Gaussian diagonal is retained.

For an active undirected graph edge $(i,j)$, v0.16 evaluates the exact cached pair and
inserts

$$
S_{ij},S_{ji},
$$

and the corresponding $2\times2$ electronic Hamiltonian blocks.

Pairs absent from the graph are set to structural zero in the sparse projected
approximation.

The matrices are stored as CSR sparse matrices.

## 5. Sparse moving-basis connection

The moving-basis matrix

$$
T_{ij}=\langle G_i|\dot G_j\rangle
$$

is evaluated for:

- all diagonal pairs;
- both orientations of every active graph edge.

It is **not** assumed Hermitian.

Metric compatibility is imposed with the sparse analogue of

$$
\boxed{
T
=
T^{(0)}
+
\frac12
\left[
\frac{S_{n+1}-S_n}{\Delta t}
-
T^{(0)}-T^{(0)\dagger}
\right].
}
$$

Thus edge appearance/disappearance enters through the discrete $\dot S$ term.

## 6. Sparse Cayley propagation

The projected equation remains

$$
iS\dot C=(H-iT)C.
$$

At the midpoint,

$$
S_m=\frac12(S_n+S_{n+1}),
$$

$$
H_m=\frac12(H_n+H_{n+1}).
$$

v0.16 solves

$$
\boxed{
\left[
S_m+\frac{\Delta t}{2}(iH_m+T_m)
\right]
C_{n+1}
=
\left[
S_m-\frac{\Delta t}{2}(iH_m+T_m)
\right]
C_n
}
$$

with `scipy.sparse.linalg.spsolve`.

No universal power-law complexity is claimed for the sparse direct solve because
fill-in depends on graph topology and ordering.

## 7. What is and is not rigorous about the sparsification

The position-only screen is a rigorous upper bound on the **overlap magnitude**.

However, a small overlap does not by itself give a universal rigorous bound on every
Hamiltonian or moving-basis matrix element because kinetic and polynomial potential
factors can amplify a small overlap.

Therefore v0.16 does **not** claim that

$$
|S_{ij}|<\tau
$$

implies a universal operator-norm error bound.

Instead, the release uses three controls:

1. conservative overlap screening;
2. the independent physical-grid TDSE defect during propagation;
3. a dense a-posteriori $S/H$ audit at the release endpoint.

## 8. Dense endpoint audit

At the final release snapshot the sparse matrices are compared against a complete
dense v0.15-style pair build.

Measured errors are

$$
\boxed{
\frac{\|S_{\rm sparse}-S_{\rm dense}\|_F}
{\|S_{\rm dense}\|_F}
=
0.0051917427
}
$$

and

$$
\boxed{
\frac{\|H_{\rm sparse}-H_{\rm dense}\|_F}
{\|H_{\rm dense}\|_F}
=
0.0039626323.
}
$$

Only

```text
3
```

off-diagonal Gaussian pairs are absent at that endpoint.

The maximum omitted exact overlap is

```text
0.02057350476995086
```

and the maximum omitted Hamiltonian $2\times2$ block Frobenius norm is

```text
0.2318916964441307.
```

The latter is why the dense audit is retained rather than inferring Hamiltonian error
from the overlap cutoff alone.

## 9. Release physics

The v0.16 representation-consistent projected-state error is

$$
\boxed{
0.0001336146
}.
$$

The original-target reduced-density error is

$$
\boxed{
0.033339541.
}
$$

The final v0.16 density matrix differs from the v0.15 dense result by

$$
\boxed{
0.00016312659
}
$$

in Frobenius norm.

Thus the sparse approximation is visible but small.

## 10. Why the compact CI benchmark is only mildly sparse

The 10–11 TBF residual-selected basis is highly overlapping by construction.

Its time-averaged off-diagonal graph sparsity is only

$$
\boxed{
5.38\%.
}
$$

Accordingly, propagation pair factorizations fall by only

$$
\boxed{
4.48\%
}
$$

relative to v0.15.

This is not presented as a large sparse speedup.

Forcing a more aggressive cutoff was tested during development and degraded the
physical benchmark.

## 11. Bounded-locality scaling benchmark

A separate chain benchmark tests the regime where Gaussian locality actually exists.

At $N=80$:

```text
active off-diagonal edges: 157
all possible off-diagonal pairs: 3160
edge fraction: 0.049684
pair-factorization reduction: 92.69 %
dense matrix assembly: 0.523068 s
sparse matrix assembly: 0.016066 s
diagnostic assembly speedup: 32.56 x
```

The fitted active-edge exponent over $N=20,40,80$ is

$$
\boxed{
1.0426
}
$$

while dense canonical pair count scales with fitted exponent

$$
\boxed{
1.9738.
}
$$

The KD-tree spatial-candidate exponent is

$$
\boxed{
1.0426.
}
$$

These exponents apply only to the synthetic bounded-locality chain.

They are not a universal claim that every Gaussian basis becomes linear-scaling.

## 12. Local-degree cost-aware adaptation

v0.15 estimated one-TBF pair cost using the global basis size $N$.

v0.16 instead predicts the candidate's local degree

$$
d_c
=
\#\left\{
i:
\text{overlap upper bound}(c,i)
\ge\tau_{\rm enter}
\right\}.
$$

Over a horizon of $h$ steps, the estimated added canonical endpoint+midpoint pair work
is

$$
\boxed{
\Delta F_c
=
2h(d_c+1).
}
$$

Thus a TBF inserted into a sparse region is cheaper than one that couples to nearly
every existing TBF.

For the release event,

```text
predicted local degree: 9
utility: 0.51124608
normalized incremental cost: 0.29137463
```

## 13. Electronic-structure cost term

v0.16 adds an optional provider-aware term to the adaptation denominator.

The geometry-cache reference model assigns different cost units to:

```text
candidate near an already cached electronic geometry
candidate requiring a new electronic geometry
```

In the deterministic demonstration:

```text
cached geometry normalized cost:
1.675

new geometry normalized cost:
3.625
```

The analytic LVC release itself has zero electronic-structure cost, so this term does
not alter its physical candidate choice.

## 14. Scientific label

The appropriate description is:

> **persistent-locality, sparse-matrix, TDSE-defect-controlled spinor-complete
> Gaussian dynamics with local-degree and provider-cost-aware basis adaptation on an
> analytic LVC benchmark.**

It is not yet:

- a production sparse AIMS implementation;
- a rigorous universal sparse-Hamiltonian truncation theorem;
- a full PySCF residual-adaptive molecular dynamics engine.
