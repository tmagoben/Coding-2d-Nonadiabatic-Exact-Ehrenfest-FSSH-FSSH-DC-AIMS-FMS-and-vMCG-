# v0.20 Detailed Derivations

## 1. Safe nuclear-overlap radius

The conservative bound is

$$
|S_{ij}^{\rm nuc}|
\le
\exp\left[-\frac12 h_{ij}r_{ij}^2\right],
\qquad
h_{ij}=\frac1{1/a_i+1/a_j}.
$$

Since $h_{ij}\ge a_{\min}/2$,

$$
|S_{ij}^{\rm nuc}|
\le
\exp\left[-\frac{a_{\min}}4r_{ij}^2\right].
$$

For geometric floor $\tau_s$, every pair outside

$$
\boxed{
R_s=\sqrt{-\frac{4\ln\tau_s}{a_{\min}}}
}
$$

is guaranteed to have nuclear overlap below $\tau_s$.

This is not a bound on the complete molecular Hamiltonian score.

## 2. Gauge-covariant centroid transport

Let

$$
O_{ci}=\Phi_c^\dagger\Phi_i.
$$

The polar link

$$
U_{ci}=\operatorname{polar}(O_{ci})
$$

transforms under local gauges as

$$
U_{ci}\rightarrow G_c^\dagger U_{ci}G_i.
$$

The local coefficient vector transforms as

$$
e_i\rightarrow G_i^\dagger e_i.
$$

Therefore the centroid vector

$$
v_i^c=U_{ci}e_i
$$

transforms as

$$
v_i^c\rightarrow G_c^\dagger v_i^c.
$$

Consequently

$$
(v_i^c)^\dagger v_j^c
$$

and

$$
(v_i^c)^\dagger H_c v_j^c
$$

are gauge invariant when $H_c\rightarrow G_c^\dagger H_cG_c$.

## 3. Pair Hermiticity

A symmetric pair centroid is used for both orientations. Thus

$$
S_{ji}=S_{ij}^*,
\qquad
H_{ji}=H_{ij}^*.
$$

v0.20 computes one off-diagonal $S/H$ block and inserts its conjugate partner.

The moving-basis seed is ordered, however:

$$
T_{ij}^{(0)}=\langle G_i|\dot G_j\rangle,
$$

so $T_{ij}^{(0)}$ and $T_{ji}^{(0)}$ are evaluated separately.

## 4. Dimensionless molecular edge score

The Hamiltonian scale is

$$
E_{ij}
=
\max\left[
\sqrt{|H_{ii}H_{jj}|},
E_{\rm floor}
\right].
$$

Then

$$
h_{ij}=\frac{|H_{ij}|}{E_{ij}},
$$

and

$$
t_{ij}
=
\Delta t
\sqrt{|T_{ij}^{(0)}|^2+|T_{ji}^{(0)}|^2}.
$$

The score is

$$
\boxed{
\eta_{ij}^2
=
(w_S|S_{ij}|)^2+
(w_Hh_{ij})^2+
(w_Tt_{ij})^2.
}
$$

It is an importance heuristic, not a rigorous operator-norm error bound.

## 5. Local omission promotion

For omitted scored candidates $D$,

$$
B_D^2=\sum_{e\in D}\eta_e^2.
$$

If $B_D>B_{\max}$, sort by decreasing $\eta_e$ and restore the largest edges until

$$
\sqrt{
B_D^2-\sum_{k=1}^r\eta_k^2
}
\le B_{\max}.
$$

At zero budget all geometrically scored candidates are restored.

## 6. Why zero budget can retain tiny dense error

The local score budget does not operate on pairs screened before electronic scoring.

Therefore $B_{\max}=0$ alone is not the full dense limit.

The formal dense limit requires

$$
B_{\max}\rightarrow0
\quad\text{and}\quad
\tau_s\rightarrow0.
$$

The release zero-budget residuals are of order $10^{-6}$ because the geometric floor
remains finite.

## 7. Metric-compatible sparse connection

The moving-basis identity is

$$
\dot S=T+T^\dagger.
$$

For a physical seed $T_0$,

$$
\boxed{
T=T_0+\frac12(\dot S-T_0-T_0^\dagger).
}
$$

Using

$$
\dot S\approx\frac{S_{n+1}-S_n}{\Delta t}
$$

makes the correction valid even when sparse edge topology changes between endpoints,
because both sparse matrices live in the same full basis dimension.

## 8. Sparse midpoint/Cayley solve

With

$$
iS\dot C=(H-iT)C
$$

and $K=iH+T$, the midpoint step is

$$
\boxed{
\left(S_m+\frac{\Delta t}{2}K_m\right)C_{n+1}
=
\left(S_m-\frac{\Delta t}{2}K_m\right)C_n.
}
$$

The implementation solves this sparse system directly and never forms $S^{-1}$.

## 9. Geometric-search controller

If a sampled omitted pair has

$$
\eta_{ij}>\eta_{\rm enter},
$$

the score threshold would retain it. Therefore the failure belongs to candidate
generation.

v0.20 changes only

$$
\tau_s\leftarrow r\tau_s,
$$

rebuilds, and re-audits.

This keeps the corrective action tied to the approximation layer that actually failed.

## 10. Dense-sentinel cache isolation

Let $\mathcal C_s$ be the production sparse electronic cache and $\mathcal C_d$ an
independent dense-validation cache.

Using the same cache would precompute all dense centroids and bias subsequent sparse
cost measurements.

The release therefore keeps

$$
\boxed{
\mathcal C_s\cap\mathcal C_d=\varnothing
}
$$

at the software cache level.

## 11. Buffered KD-tree lookup

Let $N_i$ trusted points live in the immutable KD-tree and $B$ recent insertions remain
in a direct-search buffer.

An exact query costs approximately

$$
O(\log N_i+B n_q)
$$

between rebuilds.

Rebuilding the tree costs approximately

$$
O(N_c\log N_c).
$$

The rebuild cost is explicitly counted and is not claimed to disappear.
