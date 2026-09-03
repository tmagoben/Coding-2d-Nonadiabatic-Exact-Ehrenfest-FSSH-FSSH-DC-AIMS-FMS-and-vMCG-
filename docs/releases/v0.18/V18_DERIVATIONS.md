# v0.18 Detailed Derivations

## 1. Phase alignment

For normalized wavefunctions $\psi$ and $\phi$, minimize

$$
f(\theta)
=
\|\psi-e^{i\theta}\phi\|^2.
$$

Expanding,

$$
f(\theta)
=
2
-
2\Re
\left[
e^{i\theta}
\langle\psi|\phi\rangle
\right].
$$

Write

$$
\langle\psi|\phi\rangle
=
|z|e^{i\alpha}.
$$

The maximum real overlap occurs for

$$
\theta=-\alpha.
$$

Therefore

$$
\boxed{
\phi_{\mathrm{aligned}}
=
e^{-i\arg\langle\psi|\phi\rangle}\phi.
}
$$

The minimum distance satisfies

$$
\boxed{
\epsilon_\Psi^2
=
2-2|\langle\psi|\phi\rangle|
=
2-2\sqrt F.
}
$$

This is exactly the phase convention used by v0.18.

## 2. Exact projected-reference overlap

Let

$$
|\psi_t(0)\rangle
$$

and

$$
|\psi_p(0)\rangle
$$

be target and projected initial states.

Under the same unitary propagator $U(t)$,

$$
|\psi_t(t)\rangle=U(t)|\psi_t(0)\rangle,
$$

$$
|\psi_p(t)\rangle=U(t)|\psi_p(0)\rangle.
$$

Then

$$
\begin{aligned}
\langle\psi_t(t)|\psi_p(t)\rangle
&=
\langle\psi_t(0)|
U^\dagger(t)U(t)
|\psi_p(0)\rangle
\\
&=
\boxed{
\langle\psi_t(0)|\psi_p(0)\rangle.
}
\end{aligned}
$$

Thus the exact projected-target fidelity is invariant.

This gives an unusually clean separation between initial representation error and
subsequent approximate Gaussian propagation error.

## 3. Nuclear density

For a spinor

$$
\Psi(R)
=
\sum_a\Psi_a(R)|a\rangle,
$$

tracing over the electronic state gives

$$
\boxed{
n(R)
=
\sum_a|\Psi_a(R)|^2.
}
$$

Normalization requires

$$
\int n(R)dR=1.
$$

The density L2 error is

$$
\epsilon_n
=
\left[
\int
(n_G-n_{\mathrm{ref}})^2dR
\right]^{1/2}.
$$

The total-variation distance is

$$
D_{\mathrm{TV}}
=
\frac12
\int|n_G-n_{\mathrm{ref}}|dR.
$$

## 4. Spatial moments

The centroid is

$$
\mu_\alpha
=
\int R_\alpha n(R)dR.
$$

The covariance is

$$
\Sigma_{\alpha\beta}
=
\int
(R_\alpha-\mu_\alpha)
(R_\beta-\mu_\beta)
n(R)dR.
$$

v0.18 reports

$$
\|\mu_G-\mu_{\mathrm{ref}}\|_2
$$

and

$$
\|\Sigma_G-\Sigma_{\mathrm{ref}}\|_F.
$$

These diagnostics distinguish a globally misplaced packet from an incorrectly broadened
or correlated packet.

## 5. Successive self-convergence

Assume an asymptotic discretization error

$$
u_h
=
u
+
Ch^p
+
O(h^{p+1}).
$$

Then

$$
u_h-u_{h/2}
=
Ch^p(1-2^{-p})
+
O(h^{p+1}),
$$

and

$$
u_{h/2}-u_{h/4}
=
C(h/2)^p(1-2^{-p})
+
O(h^{p+1}).
$$

Therefore

$$
\frac{
\|u_h-u_{h/2}\|
}{
\|u_{h/2}-u_{h/4}\|
}
\rightarrow
2^p.
$$

Hence

$$
\boxed{
p_{\mathrm{obs}}
=
\frac{
\ln
\left(
\|u_h-u_{h/2}\|
/
\|u_{h/2}-u_{h/4}\|
\right)
}{
\ln2
}.
}
$$

v0.18 applies this directly to phase-aligned complex Gaussian wavefunctions.

## 6. Why adaptive cadence must use physical time

Suppose enrichment is checked every $m$ numerical steps.

The physical interval is

$$
\Delta t_{\mathrm{control}}
=
m\Delta t.
$$

If $m$ is held fixed while $\Delta t$ changes, the algorithm itself changes during a
timestep convergence study.

v0.18 instead specifies

$$
\tau_{\mathrm{control}}
$$

and resolves

$$
\boxed{
m
=
\operatorname{round}
\left(
\tau_{\mathrm{control}}/\Delta t
\right).
}
$$

The same normalization is used for:

- TDSE-defect checks;
- minimum adaptation separation;
- prune age;
- sampled-audit interval;
- cost-model horizon.

For the release timestep ladder, enrichment occurs at the same physical times

$$
0.05,\ 0.10,\ 0.15
$$

for every $\Delta t$.

## 7. Batched residual ranking

Let

$$
Q
\in
\mathbb{C}^{K\times G}
$$

contain all candidate Gaussian grid values.

Earlier ranking materialized $Q$ completely.

v0.18 partitions the candidates into batches

$$
Q^{(b)}
\in
\mathbb{C}^{B\times G}.
$$

The basis-grid matrix

$$
B
\in
\mathbb{C}^{N\times G}
$$

is retained.

For one batch,

$$
X
=
B^\dagger Q^{(b)}\Delta V,
$$

and the orthogonalized candidate norm follows from the same projection algebra as the
dense implementation.

Thus arithmetic still scales over all $K$ candidates, but peak candidate storage
changes from

$$
KG
$$

to

$$
BG.
$$

The automated tests verify that batched and dense candidate rankings return the same
candidate indices and capture fractions to numerical tolerance.

## 8. Sampled omitted-edge audit

Let $D$ be the omitted edge set.

v0.18 selects

$$
J=J_{\mathrm{priority}}+J_{\mathrm{random}}
$$

audit candidates.

Priority edges are chosen near the geometric search boundary.

Random edges are selected deterministically from the omitted pair index space.

Each sampled pair receives the exact v0.17 importance score

$$
\eta_{ij}.
$$

A violation is recorded when

$$
\boxed{
\eta_{ij}
>
\eta_{\mathrm{enter}}.
}
$$

This checks whether a pair excluded by sparse search/selection would independently
qualify as an active edge.

The audit does not establish a strict bound over the unsampled remainder.

## 9. Dense sentinels

The complete dense audit remains

$$
\epsilon_H
=
\frac{
\|H_s-H_d\|_F
}{
\|H_d\|_F
}.
$$

Unlike v0.17, it is evaluated only at initial and final release sentinels.

The normal trajectory therefore avoids repeated $O(N^2)$ dense pair rebuilds while
retaining an exact beginning/end validation.

## 10. Basis-completeness interpretation

The basis ladder does not prove convergence to zero error because the maximum release
basis is still only 13 Gaussians.

It establishes the weaker but important statement

$$
\epsilon_\Psi(N=13)
<
\epsilon_\Psi(N=12)
<
\epsilon_\Psi(N=11)
<
\epsilon_\Psi(N=10)
$$

for the controlled benchmark.

That is evidence of **systematic improvement**, not proof of complete basis
convergence.
