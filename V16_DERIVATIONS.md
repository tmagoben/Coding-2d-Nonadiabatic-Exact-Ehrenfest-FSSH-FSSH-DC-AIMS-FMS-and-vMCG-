# v0.16 Detailed Derivations

## 1. Position-overlap bound

For two normalized Gaussians with real SPD width matrices $A_i$ and $A_j$, the
zero-momentum spatial overlap contains

$$
\exp\left[
-\frac12\Delta q^T
(A_i^{-1}+A_j^{-1})^{-1}
\Delta q
\right].
$$

The determinant prefactor is no larger than one.

Let

$$
a_i=\lambda_{\min}(A_i),
\qquad
a_j=\lambda_{\min}(A_j).
$$

Because

$$
A_i^{-1}
\preceq
\frac{1}{a_i}I,
\qquad
A_j^{-1}
\preceq
\frac{1}{a_j}I,
$$

we have

$$
A_i^{-1}+A_j^{-1}
\preceq
\left(
\frac1{a_i}+\frac1{a_j}
\right)I.
$$

Inverting reverses the Loewner order:

$$
(A_i^{-1}+A_j^{-1})^{-1}
\succeq
\frac{1}{
1/a_i+1/a_j
}I.
$$

Therefore

$$
\Delta q^T
(A_i^{-1}+A_j^{-1})^{-1}
\Delta q
\ge
h_{ij}\|\Delta q\|^2,
$$

with

$$
h_{ij}
=
\frac{1}{1/a_i+1/a_j}.
$$

Hence

$$
\boxed{
|S_{ij}|
\le
\exp\left[
-\frac12h_{ij}\|\Delta q\|^2
\right].
}
$$

The exact finite-momentum overlap contains an additional nonpositive real contribution,
so the same bound remains valid.

## 2. Global KD-tree radius

Let

$$
a_{\min}=\min_i a_i.
$$

Since $a_i,a_j\ge a_{\min}$,

$$
h_{ij}
=
\frac{a_ia_j}{a_i+a_j}
\ge
\frac{a_{\min}}{2}.
$$

Therefore

$$
|S_{ij}|
\le
\exp\left[
-\frac{a_{\min}}{4}\|\Delta q\|^2
\right].
$$

Requiring the right side to be at least a threshold $\tau$ gives

$$
R_\tau
=
\boxed{
\sqrt{
-\frac{4\ln\tau}{a_{\min}}
}
}.
$$

Pairs farther apart than $R_\tau$ cannot survive even the conservative overlap screen.

v0.16 uses the smaller exit threshold to define the global radius because an existing
edge is allowed to persist below the enter threshold.

## 3. Hysteretic edge rule

Let

$$
\tau_{\mathrm{exit}}<\tau_{\mathrm{enter}}.
$$

For an inactive pair,

$$
(i,j)\notin E_n,
$$

the next graph contains the edge only if

$$
|S_{ij}^{n+1}|
\ge
\tau_{\mathrm{enter}}.
$$

For an active pair,

$$
(i,j)\in E_n,
$$

the edge remains if

$$
|S_{ij}^{n+1}|
\ge
\tau_{\mathrm{exit}}.
$$

Thus the graph has memory and avoids rapid on/off switching near one cutoff.

## 4. Sparse block structure

For two electronic states per Gaussian, each active nuclear edge contributes a
$2\times2$ electronic block.

If the graph has $N$ vertices and $E$ active off-diagonal edges, then the maximum block
count is

$$
N+2E
$$

because each undirected off-diagonal edge contributes both matrix orientations.

For a dense electronic block, the nominal Hamiltonian nonzero count scales as

$$
\boxed{
4(N+2E).
}
$$

The overlap block is proportional to the $2\times2$ identity, so its actual scalar
nonzero count is smaller.

## 5. Sparse metric-compatible connection

Start from an ordered sparse seed $T^{(0)}$ and finite-difference metric derivative

$$
\dot S
\approx
\frac{S_{n+1}-S_n}{\Delta t}.
$$

Define

$$
T
=
T^{(0)}
+
\frac12
\left[
\dot S
-
T^{(0)}
-
T^{(0)\dagger}
\right].
$$

Then

$$
\boxed{
T+T^\dagger=\dot S.
}
$$

This identity remains exact in sparse matrix algebra.

## 6. Sparse Cayley equation

With

$$
K_m=iH_m+T_m,
$$

the implicit midpoint discretization is

$$
\left(
S_m+\frac{\Delta t}{2}K_m
\right)C_{n+1}
=
\left(
S_m-\frac{\Delta t}{2}K_m
\right)C_n.
$$

v0.16 stores the left matrix in sparse CSC form and solves it with `spsolve`.

The exact complexity of sparse LU depends on elimination fill-in and therefore cannot
be summarized by a single graph-independent exponent.

## 7. Local candidate degree

For a candidate Gaussian $c$, define its cheap predicted local degree

$$
d_c
=
\sum_{i=1}^{N}
\mathbf 1
\left[
B_{ci}\ge\tau_{\mathrm{enter}}
\right],
$$

where $B_{ci}$ is the conservative overlap upper bound.

This does not require exact pair factorizations.

The estimated added canonical pair work over $h$ future steps is

$$
\boxed{
\Delta F_c
=
2h(d_c+1).
}
$$

## 8. Sparse block-growth proxy

The added spinor-complete sparse block structure contains one diagonal $s\times s$
block and two oriented $s\times s$ blocks for each local neighbor.

Thus the estimated scalar nonzero growth is

$$
\boxed{
\Delta n_{\mathrm{nz}}
=
s^2(1+2d_c).
}
$$

This is a matrix-storage/work proxy.

It is not a prediction of sparse-LU fill-in.

## 9. Provider cost

Let the electronic provider return a dimensionless cost estimate

$$
C_{\mathrm{el}}(c).
$$

For the geometry-cache reference model,

$$
C_{\mathrm{el}}(c)
=
\begin{cases}
C_{\mathrm{cached}}, & d(q_c,\mathcal Q_{\mathrm{cache}})\le r_{\mathrm{reuse}},\\
C_{\mathrm{new}}, & \text{otherwise}.
\end{cases}
$$

This term is deliberately separated from Gaussian matrix cost.

In a real molecular implementation, the numerical units should be calibrated from
measured SCF/CASSCF/gradient/NAC wall times.

## 10. Local sparse utility

The v0.16 cost model combines

$$
C_c
=
\left[
r_{\mathrm{pair}}
+\frac12r_{\mathrm{sparse}}
+w_{\mathrm{el}}C_{\mathrm{el}}(c)
\right]
m_\kappa,
$$

where

$$
m_\kappa
=
1+
w_\kappa
\max
\left[
\log_{10}
\left(
\frac{\kappa_c}{\kappa_0}
\right),
0
\right].
$$

The adaptation utility is

$$
\boxed{
U_c
=
\frac{f_c}{C_c},
}
$$

where $f_c$ is the TDSE-defect capture fraction.

The residual shortlist is still the first gate.

A cheap candidate with negligible physical benefit cannot be promoted.

## 11. Sparse-operator caveat

The overlap threshold rigorously bounds $|S_{ij}|$.

It does **not** produce a universal bound on

$$
\|H_{ij}\|,
\qquad
\|T_{ij}\|.
$$

For polynomial Hamiltonians these elements are overlap multiplied by coordinate and
momentum moments. Those prefactors can be larger than one.

Therefore v0.16 explicitly audits the final sparse $S$ and $H$ against dense matrices.

This is an a-posteriori validation strategy, not a sparsification theorem.

## 12. Bounded-locality graph scaling

For a chain with approximately constant local neighbor count,

$$
E=O(N).
$$

Then active pair algebra is

$$
O(Nd^3),
$$

while dense canonical pair algebra is

$$
O(N^2d^3).
$$

The KD-tree neighbor query is approximately

$$
O(N\log N+E)
$$

for bounded local density.

Worst-case behavior remains quadratic when the safe global radius includes nearly the
entire basis.

That worst-case behavior is exactly what occurs in the compact CI benchmark: the basis
is highly overlapping and is not a sparse asymptotic regime.
