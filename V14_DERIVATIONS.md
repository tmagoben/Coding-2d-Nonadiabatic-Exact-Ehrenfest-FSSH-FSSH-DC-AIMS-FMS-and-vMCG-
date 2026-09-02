# v0.14 Detailed Derivations

## 1. Time-dependent basis dimension

The approximation is

$$
|\Psi(t)\rangle
=
\sum_{\mu=1}^{M(t)}
C_\mu(t)|\Xi_\mu(t)\rangle.
$$

For the two-state spinor-complete basis,

$$
M(t)=2N(t).
$$

Between discrete adaptation events, $N(t)$ is constant and the ordinary moving-basis
derivation applies.

At an adaptation time $t_*$, the basis changes discontinuously but the represented
wavefunction is constrained to remain continuous.

---

## 2. Growth insertion continuity

Let the old coefficient vector be

$$
C^-.
$$

Adding one nuclear Gaussian with two electronic components gives

$$
\boxed{
C^+
=
\begin{pmatrix}
C^-\\
0\\
0
\end{pmatrix}.
}
$$

If the new basis functions are

$$
|\Xi_{M+1}\rangle,
\qquad
|\Xi_{M+2}\rangle,
$$

then

$$
\begin{aligned}
|\Psi^+\rangle
&=
\sum_{\mu=1}^{M}
C_\mu^-|\Xi_\mu\rangle
+
0|\Xi_{M+1}\rangle
+
0|\Xi_{M+2}\rangle
\\
&=
\boxed{
|\Psi^-\rangle.
}
\end{aligned}
$$

No renormalization is required in exact arithmetic after zero-amplitude insertion.

---

## 3. TDSE defect

The exact TDSE is

$$
i|\dot\Psi\rangle=H|\Psi\rangle.
$$

For the finite moving basis define

$$
\boxed{
|\mathcal R\rangle
=
i|\dot\Psi_G\rangle-H|\Psi_G\rangle.
}
$$

The normalized monitor is

$$
\boxed{
\eta
=
\frac{\|\mathcal R\|}{\|H\Psi_G\|}.
}
$$

The denominator makes the threshold less dependent on an arbitrary overall energy
scale.

It is still a benchmark-specific error indicator rather than a universal error bound
on every observable.

---

## 4. Candidate orthogonalization

For candidate nuclear Gaussian $g_c$ and current nuclear basis $g_i$,

$$
S_{ij}=\langle g_i|g_j\rangle,
$$

$$
s_i=\langle g_i|g_c\rangle.
$$

The projection coefficients satisfy

$$
S\alpha=s.
$$

Hence

$$
g_c^\perp
=
g_c-\sum_i\alpha_i g_i.
$$

Its squared norm is

$$
\boxed{
n_c
=
\langle g_c^\perp|g_c^\perp\rangle
=
1-s^\dagger S^{-1}s.
}
$$

---

## 5. TDSE-defect capture

Write

$$
\mathcal R
=
\sum_a\mathcal R_a|d_a\rangle.
$$

The best correction in the new spinor-complete candidate pair is

$$
\delta\Psi_c
=
g_c^\perp
\sum_a\beta_a|d_a\rangle,
$$

with

$$
\boxed{
\beta_a
=
\frac{\langle g_c^\perp|\mathcal R_a\rangle}{n_c}.
}
$$

The squared defect norm captured is

$$
\boxed{
\Delta_c^{\rm TDSE}
=
\frac{
\sum_a
|\langle g_c^\perp|\mathcal R_a\rangle|^2
}{
n_c
}.
}
$$

The capture fraction is

$$
\boxed{
f_c
=
\frac{\Delta_c^{\rm TDSE}}{\|\mathcal R\|^2}.
}
$$

The release requires

$$
f_c
$$

to exceed a configurable minimum before a candidate is accepted.

---

## 6. Hysteretic control rule

Define

$$
\eta_{\rm add}
>
\eta_{\rm remove}.
$$

The control law is

$$
\boxed{
\begin{cases}
\text{consider enrichment},
&
\eta\ge\eta_{\rm add},
\\[4pt]
\text{no basis-size action},
&
\eta_{\rm remove}<\eta<\eta_{\rm add},
\\[4pt]
\text{accumulate pruning patience},
&
\eta\le\eta_{\rm remove}.
\end{cases}
}
$$

This is not quantum-mechanical physics.

It is numerical control logic wrapped around the quantum propagation.

That distinction is explicit in the implementation.

---

## 7. Exact leave-one-out pruning loss

Let the nuclear basis be partitioned into candidate $j$ and all remaining functions.

Write

$$
g_j
=
P_{-j}g_j
+
h_j,
$$

where

$$
h_j
=
(I-P_{-j})g_j.
$$

By construction,

$$
h_j\perp\operatorname{span}\{g_i:i\neq j\}.
$$

The current spinor-complete wavefunction is

$$
\Psi
=
\sum_{i\neq j}g_i\mathbf C_i
+
g_j\mathbf C_j.
$$

Substitute the decomposition:

$$
\Psi
=
\left[
\sum_{i\neq j}g_i\mathbf C_i
+
(P_{-j}g_j)\mathbf C_j
\right]
+
h_j\mathbf C_j.
$$

The bracketed term belongs to the retained basis.

The final term is exactly orthogonal to it.

Therefore the squared best-projection loss after deleting $j$ is

$$
\boxed{
L_j
=
\|h_j\|^2
\|\mathbf C_j\|^2.
}
$$

---

## 8. Schur-complement identity

Let

$$
S
$$

be the full nuclear Gram matrix.

For a positive-definite Gram matrix,

$$
\boxed{
\|h_j\|^2
=
\frac{1}{(S^{-1})_{jj}}.
}
$$

This follows from the Schur complement of the retained block.

Hence

$$
\boxed{
L_j
=
\frac{
\displaystyle\sum_a|C_{ja}|^2
}{
(S^{-1})_{jj}
}.
}
$$

The fractional pruning loss is

$$
\boxed{
\ell_j
=
\frac{L_j}{\langle\Psi|\Psi\rangle}.
}
$$

All $L_j$ are obtained from one inverse/solve.

---

## 9. Why this pruning score is better than coefficient magnitude alone

A small coefficient does not necessarily imply a removable Gaussian in a
nonorthogonal basis.

Likewise, a moderately large coefficient can multiply a basis direction that is
almost completely represented by its neighbors.

The factor

$$
\frac{1}{(S^{-1})_{jj}}
$$

measures the genuinely independent norm of basis direction $j$.

Thus the deletion score combines:

1. coefficient importance;
2. nonorthogonal geometric redundancy.

---

## 10. Hermitian half-build derivation

For Hermitian operators,

$$
S_{ji}
=
\langle g_j|g_i\rangle
=
S_{ij}^*,
$$

and

$$
H_{ji}
=
\langle\Xi_j|H|\Xi_i\rangle
=
H_{ij}^\dagger.
$$

Therefore only pairs

$$
i\le j
$$

need direct evaluation.

The number is

$$
\sum_{i=1}^N(N-i+1)
=
\boxed{
\frac{N(N+1)}{2}.
}
$$

The old ordered-pair count is

$$
N^2.
$$

The relative reduction is

$$
\boxed{
R_N
=
1-\frac{N(N+1)/2}{N^2}
=
\frac12-\frac{1}{2N}.
}
$$

As

$$
N\to\infty,
$$

$$
R_N\to\frac12.
$$

---

## 11. Why the moving-basis T matrix cannot be half-built the same way

The matrix

$$
T_{ij}
=
\langle\Xi_i|\dot\Xi_j\rangle
$$

is not Hermitian.

Instead,

$$
\boxed{
\dot S=T+T^\dagger.
}
$$

Therefore

$$
T_{ji}\neq T_{ij}^*
$$

in general.

v0.14 does not incorrectly impose Hermiticity on $T$.

This is why the moving-basis $T$ construction remains an

$$
O(N^2)
$$

ordered-pair operation and is a major observed runtime cost.

---

## 12. Dense Cayley solve complexity

The spinor-complete electronic matrix dimension is

$$
M=sN.
$$

The Cayley step solves a dense system

$$
A C_{n+1}=b.
$$

A generic dense factorization requires

$$
\boxed{
O(M^3)
=
O(s^3N^3)
}
$$

time and

$$
\boxed{
O(M^2)
=
O(s^2N^2)
}
$$

memory.

At the release peak,

$$
s=2,
\qquad
N=11,
\qquad
M=22.
$$

Thus the cubic solve is asymptotically important but numerically tiny at the present
benchmark size.

---

## 13. Candidate-ranking matrix form

Let

$$
B\in\mathbb C^{N\times G}
$$

contain current nuclear Gaussian values on the diagnostic grid.

Let

$$
Q\in\mathbb C^{K\times G}
$$

contain candidate Gaussians.

Then

$$
S
=
BB^\dagger\Delta A,
$$

and

$$
X
=
BQ^\dagger\Delta A
$$

contains candidate overlaps with the current span.

Solve

$$
\boxed{
A=S^{-1}X.
}
$$

The orthogonal norms are

$$
\boxed{
n_k
=
\|q_k\|^2
-
x_k^\dagger a_k.
}
$$

---

## 14. Residual contractions

Let

$$
R\in\mathbb C^{G\times s}
$$

contain the TDSE residual electronic components.

Compute

$$
B_R=B^*R\Delta A,
$$

$$
Q_R=Q^*R\Delta A.
$$

Then for every candidate simultaneously,

$$
\boxed{
b_k^\perp
=
(Q_R)_k
-
a_k^\dagger B_R.
}
$$

Hence

$$
\boxed{
\Delta_k^{\rm TDSE}
=
\frac{\|b_k^\perp\|^2}{n_k}.
}
$$

This is why the prepared implementation can rank hundreds of candidates with dense
matrix products instead of Python-level nested integration loops.

---

## 15. Candidate-ranking complexity

The dominant contractions scale as

$$
O(KGN)
$$

and

$$
O(KGs).
$$

The overlap solve for all candidates contributes

$$
O(N^3+N^2K).
$$

Thus

$$
\boxed{
T_{\rm rank}
=
O\left(
KG(N+s)+N^2K+N^3
\right).
}
$$

This is one of the highest-complexity v0.14 components when $K$ and $G$ are large.

---

## 16. Defect-evaluation complexity

One defect check requires:

1. moving-basis $T$ construction:
   $O(N^2d^3)$ in the dense unequal-width implementation;
2. projected coefficient solve:
   $O((sN)^3)$;
3. $\Psi$ and $\dot\Psi$ reconstruction:
   $O(NGs)$;
4. FFT kinetic action:
   $O(sG\log G)$;
5. potential action:
   $O(Gs^2)$.

Therefore

$$
\boxed{
T_{\rm defect}
=
O\left(
N^2d^3
+
(sN)^3
+
NGs
+
sG\log G
+
Gs^2
\right).
}
$$

---

## 17. Total adaptive propagation complexity

For $T$ time steps and a defect interval $m$,

$$
Q\approx\frac{T}{m}
$$

defect checks are performed.

Ignoring basis-size variation for notation, the total reference scaling is roughly

$$
\boxed{
T_{\rm total}
=
T\left[
O(N^2d^3)
+
O((sN)^3)
\right]
+
\frac{T}{m}
\left[
T_{\rm defect}
+
T_{\rm candidate}
\right].
}
$$

Adaptation can reduce $N$ but adds defect-monitoring and candidate-search overhead.

This is the fundamental accuracy/cost tradeoff.

---

## 18. Memory complexity

The dense electronic matrices require

$$
O((sN)^2)
$$

complex values.

The dynamic candidate grid requires

$$
O(KG).
$$

Current-basis grid reconstruction requires

$$
O(NGs).
$$

Hence

$$
\boxed{
M_{\rm total}
=
O(
s^2N^2+KG+NGs
).
}
$$

In the release candidate-ranking event,

$$
K=560,
\qquad
G=1600.
$$

The candidate grid alone therefore contains

$$
896000
$$

complex numbers, or about

$$
13.7\;{\rm MiB}
$$

for `complex128`.

This is larger than the tiny $22\times22$ dense quantum matrices in the current
benchmark.

---

## 19. Current-runtime versus asymptotic complexity

The release timing shows that the pairwise Gaussian algebra dominates at $N\approx11$.

That does not contradict the

$$
O((sN)^3)
$$

dense-solve scaling.

The matrix dimension is only 22, so the cubic solve has a very small absolute constant
here.

As $N$ grows, the dense solve increases cubically while pairwise Gaussian algebra grows
quadratically in $N$.

Thus the dominant bottleneck can change with problem size.

---

## 20. Complexity of pruning

All one-Gaussian deletion scores require one solve/inverse of the nuclear overlap
matrix:

$$
O(N^3).
$$

The coefficient-row scoring is then

$$
O(Ns).
$$

If a Gaussian is actually removed, projection into the retained basis requires an
additional dense solve for each electronic component.

Hence a pruning audit is

$$
\boxed{
O(N^3+Ns),
}
$$

with an accepted pruning event still dominated by

$$
O(N^3).
$$
