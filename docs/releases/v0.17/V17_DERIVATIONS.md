# v0.17 Detailed Derivations

## 1. Hamiltonian-normalized edge term

The off-diagonal block $H_{ij}$ has energy units.

To combine it with the dimensionless overlap, v0.17 defines

$$
E_{ij}
=
\max
\left[
\sqrt{
\|H_{ii}\|_F
\|H_{jj}\|_F
},
E_{\mathrm{floor}}
\right].
$$

Then

$$
\boxed{
h_{ij}
=
\frac{\|H_{ij}\|_F}{E_{ij}}
}
$$

is dimensionless.

The geometric-mean diagonal scale is symmetric under $i\leftrightarrow j$.

## 2. Time-connection term

The moving-basis matrix element has inverse-time units in atomic units.

Multiplication by the integration step gives a dimensionless quantity.

Because $T$ is ordered and generally non-Hermitian, v0.17 uses both orientations:

$$
\boxed{
t_{ij}
=
\Delta t
\sqrt{
|T_{ij}|^2+|T_{ji}|^2
}.
}
$$

No Hermitian assumption is introduced.

## 3. Combined score

The score is

$$
\eta_{ij}
=
\sqrt{
(w_S|S_{ij}|)^2
+
(w_Hh_{ij})^2
+
(w_Tt_{ij})^2
}.
$$

It is nonnegative, symmetric under consistent pair evaluation, and prevents one
channel from cancelling another.

Setting $w_H=w_T=0$ recovers an overlap-only score.

## 4. Global omitted-score proxy

Let $D$ be locally scored edges tentatively removed by the per-edge hysteresis rule.

Define

$$
B_{\mathrm{local}}^2
=
\sum_{e\in D}\eta_e^2.
$$

If $B_{\mathrm{local}}>B_{\max}$, sort omitted edges so that

$$
\eta_1\ge\eta_2\ge\cdots.
$$

Promote edges in that order until

$$
\sqrt{
B_{\mathrm{local}}^2
-
\sum_{k=1}^{m}\eta_k^2
}
\le
B_{\max}.
$$

This is exactly what the implementation performs.

## 5. Why the budget is only a proxy

In general,

$$
\left\|
\sum_e\Delta H_e
\right\|
$$

is not equal to

$$
\sqrt{
\sum_e \eta_e^2
}.
$$

Edge contributions can overlap in matrix support, sparse-LU fill can change, and the
moving nonorthogonal metric couples blocks.

Therefore $B_{\mathrm{local}}$ only prevents uncontrolled accumulation of small local
omissions.

The dense audit remains the authoritative matrix-error check.

## 6. Dense audit relaxation

Suppose at audit $k$ either

$$
\epsilon_S>\epsilon_S^{\max}
$$

or

$$
\epsilon_H>\epsilon_H^{\max}.
$$

With $0<r<1$,

$$
\eta_{\mathrm{enter}}^{(k+1)}
=
r\eta_{\mathrm{enter}}^{(k)},
$$

$$
\eta_{\mathrm{exit}}^{(k+1)}
=
r\eta_{\mathrm{exit}}^{(k)},
$$

and

$$
\tau_{\mathrm{search}}^{(k+1)}
=
r\tau_{\mathrm{search}}^{(k)}.
$$

Because all thresholds only decrease, the controller is one-sided.

## 7. Threshold-convergence expectation

For a fixed snapshot, decreasing the enter threshold restores exact dense matrix
blocks.

If the edge sets are nested,

$$
E_1\subseteq E_2\subseteq\cdots,
$$

then for Frobenius norm

$$
\|A-A_{E_{k+1}}\|_F
\le
\|A-A_{E_k}\|_F.
$$

The release explicitly tests this for both $S$ and $H$.

## 8. Budget-convergence expectation

For fixed score threshold, tightening $B_{\max}$ promotes the largest remaining
omitted scores.

At $B_{\max}=0$, every locally scored candidate edge is restored.

In the release final snapshot, the strictest budget row produces zero audited S/H error
because every dense pair is restored.

## 9. Audit amortization

Let one dense audit cost

$$
C_{\mathrm{audit}}
=
O(N^2d^3+s^2N^2).
$$

If audits are performed every $m$ propagation steps, the amortized audit work per step
is

$$
\boxed{
C_{\mathrm{audit}}^{\mathrm{amortized}}
=
\frac{1}{m}
O(N^2d^3+s^2N^2).
}
$$

v0.17 uses $m=20$.

This is explicitly a temporary correctness cost.

## 10. Construction complexity matters separately from storage sparsity

A sparse matrix can still be expensive to construct if every pair is examined.

v0.17 therefore measures:

- KD-tree spatial candidates;
- pair-specific geometric screens;
- exact S/H/T pair checks;
- active graph edges;
- dense canonical pair count.

For a bounded-locality chain, both exact scoring candidates and active edges approach
linear growth, while the dense canonical count remains quadratic.

Worst-case dense configurations remain $O(N^2)$.
