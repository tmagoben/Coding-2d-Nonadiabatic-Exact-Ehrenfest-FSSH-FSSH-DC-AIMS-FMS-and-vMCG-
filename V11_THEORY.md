# v0.11 Theory: Basis Completeness, Optimized Spawning, and Unequal-Width Gaussian Dynamics

Version 0.10 established an important negative result.

The stronger two-dimensional conical-intersection passage conserved the generalized
Gaussian norm, but the small spawned basis did **not** reproduce the exact reduced
electronic density.

That observation determines the purpose of v0.11:

$$
\boxed{
\text{improve the adaptive Gaussian basis itself}
}
$$

rather than changing the exact benchmark, weakening the acceptance criteria, or adding
another unrelated formalism.

The scientific ladder is now

$$
\boxed{
\text{v0.10 convergence diagnosis}
\rightarrow
\text{unequal-width Gaussian algebra}
\rightarrow
\text{multi-generation spawning}
\rightarrow
\text{constrained child optimization}
\rightarrow
\text{width-diverse basis enrichment}
\rightarrow
\text{basis-completeness campaign}.
}
$$

Atomic units are used throughout.

---

# 1. What failed in v0.10

For the default strong-CI passage,

$$
q_0=(-0.60,0.25),
\qquad
p_0=(10,0),
\qquad
M=5,
\qquad
t_f=0.60,
$$

v0.10 found that a small graph-Gaussian basis remained much too close to a pure
electronic state.

The exact wavepacket developed substantial electron-nuclear correlation, while the
Gaussian basis did not branch sufficiently in nuclear phase space.

The important conclusion was not merely

> the populations are wrong.

It was more specific:

$$
\boxed{
\text{the adaptive basis does not span enough of the branched wavepacket support.}
}
$$

That is a basis-completeness problem.

---

# 2. Why adding more copies of the same Gaussian is not enough

Suppose two basis functions have nearly identical

$$
(q,p,A)
$$

and the same electronic state.

Then

$$
|\langle g_i|g_j\rangle|
\approx1.
$$

The overlap matrix develops a small eigenvalue and

$$
\kappa(S)
\gg1.
$$

Adding such a function increases the nominal TBF count without adding a genuinely new
direction to the Hilbert-space span.

Therefore v0.11 distinguishes

$$
\boxed{
N_{\rm TBF}
}
$$

from

$$
\boxed{
\text{effective independent basis dimension}.
}
$$

A useful new TBF must add phase-space or shape information that is not already
represented.

---

# 3. Standard frozen-Gaussian AIMS versus the v0.11 extension

AIMS and FMS are conventionally formulated with multidimensional frozen Gaussian
trajectory basis functions.

In a standard implementation, widths are generally fixed during propagation and are
often common to the TBFs.

v0.11 preserves the **frozen-within-each-TBF** idea:

$$
\boxed{
\dot A_i=0
}
$$

during ordinary propagation.

However, different spawned TBFs are now allowed to receive different real
positive-definite width matrices.

Thus

$$
\boxed{
A_i\neq A_j
}
$$

is permitted.

This is a deliberate basis-completeness extension.

It should not be described as standard AIMS width propagation.

It is also not vMCG: in vMCG-type methods, Gaussian parameters can evolve
variationally as part of the equations of motion.

v0.11 instead uses a **discrete width bank at basis creation** while keeping each
created Gaussian frozen afterward.

---

# 4. General unequal-width Gaussian

For TBF $i$,

$$
\boxed{
g_i(q)
=
N_i
\exp
\left[
-\frac12(q-q_i)^TA_i(q-q_i)
+
ip_i^T(q-q_i)
\right],
}
$$

where

$$
A_i=A_i^T>0.
$$

The normalization is

$$
\boxed{
N_i
=
\left(
\frac{\det A_i}{\pi^D}
\right)^{1/4}.
}
$$

Unlike the earlier equal-width formulas, v0.11 permits

$$
A_i\neq A_j.
$$

---

# 5. Product of two unequal-width Gaussians

The cross density is

$$
g_i^*(q)g_j(q).
$$

Define

$$
\boxed{
B_{ij}=A_i+A_j.
}
$$

The linear coefficient in the exponent is

$$
\boxed{
\ell_{ij}
=
A_iq_i+A_jq_j+i(p_j-p_i).
}
$$

The complex cross centroid is therefore

$$
\boxed{
\mu_{ij}
=
B_{ij}^{-1}\ell_{ij}.
}
$$

The cross covariance is

$$
\boxed{
\Sigma_{ij}
=
B_{ij}^{-1}.
}
$$

These two objects replace the equal-width expressions

$$
\frac{q_i+q_j}{2}
+
\frac{i}{2}A^{-1}(p_j-p_i)
$$

and

$$
\frac12A^{-1}.
$$

---

# 6. Exact unequal-width overlap

The exact normalized overlap is

$$
S_{ij}
=
\langle g_i|g_j\rangle.
$$

After completing the square,

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
\exp
\left[
c_{ij}
+
\frac12
\ell_{ij}^T
B_{ij}^{-1}
\ell_{ij}
\right],
}
$$

where

$$
c_{ij}
=
-\frac12q_i^TA_iq_i
-\frac12q_j^TA_jq_j
+
ip_i^Tq_i
-
ip_j^Tq_j.
$$

If

$$
A_i=A_j=A,
$$

the prefactor reduces to one and the result collapses to the equal-width formula used
in earlier releases.

That reduction is unit tested.

---

# 7. Real overlap saddle point

The maximum of

$$
|g_i(q)g_j(q)|
$$

is independent of momentum phase.

The real saddle point satisfies

$$
\nabla_q
\left[
\frac12(q-q_i)^TA_i(q-q_i)
+
\frac12(q-q_j)^TA_j(q-q_j)
\right]
=0.
$$

Therefore

$$
(A_i+A_j)q_s
=
A_iq_i+A_jq_j.
$$

Hence

$$
\boxed{
q_s
=
(A_i+A_j)^{-1}
(A_iq_i+A_jq_j).
}
$$

For equal widths,

$$
q_s=\frac{q_i+q_j}{2}.
$$

This weighted saddle is required for a mathematically consistent SPA0/SPA1 treatment
of width-diverse TBF pairs.

---

# 8. Unequal-width gradient matrix element

For Gaussian $j$,

$$
\nabla g_j
=
[-A_j(q-q_j)+ip_j]g_j.
$$

Since the first cross moment is

$$
\langle q\rangle_{ij}
=
\mu_{ij},
$$

we obtain

$$
\boxed{
G_{ij}
\equiv
\langle g_i|\nabla g_j\rangle
=
[-A_j(\mu_{ij}-q_j)+ip_j]S_{ij}.
}
$$

This formula is exact for real positive-definite $A_i$ and $A_j$.

---

# 9. Unequal-width kinetic matrix element

Let

$$
M
$$

be the generalized-coordinate mass matrix and

$$
B_M=M^{-1}.
$$

The kinetic operator is

$$
\hat T
=
-\frac12\nabla^TB_M\nabla.
$$

Using integration by parts,

$$
T_{ij}
=
\frac12
\langle
\nabla g_i
|
B_M
|
\nabla g_j
\rangle.
$$

Define

$$
u_i
=
-A_i(\mu_{ij}-q_i)-ip_i,
$$

$$
u_j
=
-A_j(\mu_{ij}-q_j)+ip_j.
$$

The fluctuation contribution is

$$
\boxed{
\operatorname{Tr}
\left[
A_iB_MA_j\Sigma_{ij}
\right].
}
$$

Therefore

$$
\boxed{
T_{ij}
=
\frac12S_{ij}
\left[
u_i^TB_Mu_j
+
\operatorname{Tr}
(
A_iB_MA_j\Sigma_{ij}
)
\right].
}
$$

When

$$
A_i=A_j=A,
$$

$$
\Sigma_{ij}=\frac12A^{-1},
$$

and the fluctuation term becomes

$$
\frac12\operatorname{Tr}(B_MA),
$$

recovering the earlier expression.

---

# 10. Moving-basis connection with unequal widths

The v0.11 TBF widths are frozen during ordinary propagation, but the algebra is written
more generally so that

$$
\dot A_j
$$

can be supplied.

For

$$
g_j
=
N_j
e^{
-\frac12\xi_j^TA_j\xi_j
+
ip_j^T\xi_j
},
$$

where

$$
\xi_j=q-q_j,
$$

differentiate the logarithm:

$$
\frac{\dot g_j}{g_j}
=
\frac14\operatorname{Tr}(A_j^{-1}\dot A_j)
+
(A_j\xi_j-ip_j)^T\dot q_j
+
i\xi_j^T\dot p_j
-
\frac12\xi_j^T\dot A_j\xi_j.
$$

Using the cross first and second moments gives

$$
\boxed{
\begin{aligned}
\langle g_i|\dot g_j\rangle
=
S_{ij}
\Big[
&
\frac14\operatorname{Tr}(A_j^{-1}\dot A_j)
\\
&+
(A_jy_{ij}-ip_j)^T\dot q_j
+
iy_{ij}^T\dot p_j
\\
&-
\frac12
\left(
y_{ij}^T\dot A_jy_{ij}
+
\operatorname{Tr}(\dot A_j\Sigma_{ij})
\right)
\Big],
\end{aligned}
}
$$

where

$$
y_{ij}=\mu_{ij}-q_j.
$$

For the present v0.11 propagator,

$$
\dot A_j=0.
$$

The more general expression is retained so a future variational-width implementation
does not require replacing the algebra again.

---

# 11. SPA0/SPA1 with unequal widths

The electronic quantity is expanded around the real overlap saddle

$$
q_s.
$$

For scalar or matrix-valued smooth field $F$,

$$
F(q)
\approx
F(q_s)
+
\nabla F(q_s)\cdot(q-q_s).
$$

The Gaussian first moment gives

$$
\int
g_i^*
F(q)
g_jdq
\approx
S_{ij}
\left[
F(q_s)
+
\nabla F(q_s)\cdot(\mu_{ij}-q_s)
\right].
$$

Thus

$$
\boxed{
\text{SPA0}:
\quad
F(q)\rightarrow F(q_s),
}
$$

and

$$
\boxed{
\text{SPA1}:
\quad
F(q)\rightarrow
F(q_s)
+
\nabla F(q_s)\cdot(\mu_{ij}-q_s).
}
$$

The old arithmetic midpoint is no longer used when the two widths differ.

---

# 12. Why child placement matters

The 2009 optimal-spawning work identified a central efficiency problem:

a spawned child is useful only if it can acquire amplitude from the parent through a
significant off-diagonal Hamiltonian coupling.

A child that satisfies energy conservation but has negligible coupling may enlarge the
basis without improving the represented wavefunction.

The conceptual optimization is therefore

$$
\boxed{
\max_{q_c,p_c}
\left|
H_{pc}
\right|
}
$$

subject to

$$
\boxed{
E_{\rm child}^{\rm class}
=
E_{\rm parent}^{\rm class}.
}
$$

The full production optimal-spawning algorithm is a nonlinear constrained
optimization.

v0.11 implements a deliberately simpler **local discrete search inspired by this
principle**.

---

# 13. Parent energy

For parent state $a$,

$$
\boxed{
\mathcal E_p
=
E_a(q_p)
+
\frac12p_p^TM^{-1}p_p.
}
$$

A candidate child on state $b$ at position $q_c$ must satisfy

$$
\boxed{
E_b(q_c)
+
\frac12p_c^TM^{-1}p_c
=
\mathcal E_p.
}
$$

This equality is enforced for every accepted spawn candidate.

---

# 14. Candidate child position

v0.11 searches a small local set of displacements

$$
q_c
=
q_p+s\hat n.
$$

The candidate direction $\hat n$ can be chosen from:

1. the nonadiabatic coupling vector;
2. the target-state force direction;
3. the parent momentum direction.

The displacement values are a symmetric local grid such as

$$
s\in
\{0,\pm0.05\}.
$$

This is a transparent finite search, not a black-box optimizer.

---

# 15. Energy-conserving momentum at a shifted child position

Choose a momentum-adjustment direction

$$
\hat m.
$$

Write

$$
p_c=p_p+\lambda\hat m.
$$

Let

$$
B=M^{-1}.
$$

The child energy constraint becomes

$$
\frac12
(p_p+\lambda\hat m)^T
B
(p_p+\lambda\hat m)
+
E_b(q_c)
=
\mathcal E_p.
$$

Expand:

$$
a\lambda^2+2b\lambda+c=0,
$$

with

$$
\boxed{
a=\hat m^TB\hat m,
}
$$

$$
\boxed{
b=p_p^TB\hat m,
}
$$

$$
\boxed{
c=p_p^TBp_p
+
2[E_b(q_c)-\mathcal E_p].
}
$$

The discriminant is

$$
\boxed{
D=b^2-ac.
}
$$

If

$$
D<0,
$$

the candidate cannot satisfy the local classical energy constraint along that momentum
direction and is discarded.

If

$$
D\ge0,
$$

$$
\boxed{
\lambda
=
\frac{-b\pm\sqrt D}{a}.
}
$$

v0.11 selects the real root with the smaller $|\lambda|$.

---

# 16. Width bank

For accepted phase-space locations, v0.11 considers

$$
\boxed{
A_c=s_AA_p
}
$$

with a small positive scale bank such as

$$
s_A\in
\{0.65,1.0,1.55\}.
$$

Interpretation:

- $s_A<1$: broader coordinate-space Gaussian;
- $s_A=1$: inherited width;
- $s_A>1$: narrower coordinate-space Gaussian.

Because the Gaussian is normalized,

$$
\operatorname{Cov}(q)
=
\frac12A^{-1}.
$$

Thus increasing $A$ narrows the packet in position and broadens it in momentum.

---

# 17. Local coupling proxy

At the common saddle frame of parent state $a$ and target state $b$, the adiabatic
electronic Hamiltonian is diagonal at zeroth order.

The leading local off-diagonal electronic Taylor term is

$$
F_{ab,\alpha}
=
(E_b-E_a)d_{ab,\alpha}.
$$

The SPA1 parent-child electronic coupling proxy is

$$
\boxed{
V_{pc}^{(1)}
=
S_{pc}^{\rm nuc}
\,
\mathbf F_{ab}(q_s)
\cdot
(\mu_{pc}-q_s).
}
$$

v0.11 ranks candidates by

$$
\boxed{
\left|
V_{pc}^{(1)}
\right|.
}
$$

This is not claimed to be the exact full off-diagonal AIMS Hamiltonian element.

It is a local, gauge-consistent first-order proxy that is cheap enough to evaluate for
many candidate children.

---

# 18. Novelty penalty

A child should also add a genuinely new basis direction.

For candidate $c$, define its maximum overlap with existing TBFs on the target state:

$$
\boxed{
S_{\max}
=
\max_{k\in b}
|\langle g_k|g_c\rangle|.
}
$$

If

$$
S_{\max}\ge S_{\rm block},
$$

the child is rejected.

For ranking, v0.11 defines a smooth novelty factor

$$
\boxed{
\nu
=
(1-S_{\max}^2)^\beta,
}
$$

with default

$$
\beta=\frac12.
$$

The ranking score is

$$
\boxed{
\mathcal J
=
|V_{pc}^{(1)}|\nu.
}
$$

This is a v0.11 heuristic extension.

It explicitly balances coupling strength against basis redundancy.

---

# 19. Multiple children per spawning event

One optimized child may still fail to represent a bifurcating nuclear wavepacket.

v0.11 can therefore select the top

$$
N_{\rm child}
$$

nonredundant candidates from the ranked local search.

Sibling candidates are required to satisfy

$$
\boxed{
|\langle g_c^{(1)}|g_c^{(2)}\rangle|
<
S_{\rm sibling}.
}
$$

Every child enters with

$$
\boxed{
C_{\rm child}=0.
}
$$

Therefore adding the basis functions does not discontinuously alter the instantaneous
represented wavefunction.

Amplitude appears only from subsequent coupled propagation.

---

# 20. Multi-generation spawning

Let generation zero be the initial TBF set.

If child $c$ is spawned from parent $p$,

$$
\boxed{
G_c=G_p+1.
}
$$

Unlike a one-generation teaching model, v0.11 permits descendants to become parents
themselves:

$$
0
\rightarrow
1
\rightarrow
2
\rightarrow
3
\rightarrow\cdots.
$$

A configurable

$$
G_{\max}
$$

prevents uncontrolled recursion.

The release stores the full lineage:

```text
uid
parent_uid
generation
birth step
birth time
birth state
width scale from parent
```

so every branch can be reconstructed after the run.

---

# 21. Why multi-generation spawning matters near a CI

A wavepacket can:

1. transfer to another electronic state;
2. propagate on the new surface;
3. encounter another strong-coupling region;
4. transfer again;
5. partially recombine with earlier branches.

A basis that permits only

$$
\text{initial parent}\rightarrow\text{one child}
$$

cannot systematically represent this hierarchy.

Multi-generation spawning creates a basis tree that can represent repeated electronic
branching.

The coefficients remain globally coupled, so the resulting state is not a classical
tree of mutually independent trajectories.

---

# 22. Basis growth versus conditioning

Larger basis size is useful only until numerical redundancy dominates.

For overlap matrix

$$
S,
$$

v0.11 monitors

$$
\boxed{
\kappa(S)
=
\frac{\lambda_{\max}}{\lambda_{\min}}.
}
$$

The inherited projection-based pruning remains available.

If a basis function is removed, the retained coefficient vector solves

$$
\boxed{
S_{KK}C'
=
S_{K,\mathrm{all}}C.
}
$$

The wavefunction projection loss is measured explicitly.

v0.11 does not silently renormalize away pruning loss.

---

# 23. Canonical effective occupation of the nonorthogonal basis

Diagonalize

$$
S=U\Lambda U^\dagger.
$$

The canonical orthonormal basis is

$$
\boxed{
\chi
=
\Phi U\Lambda^{-1/2}.
}
$$

If

$$
\Psi=\Phi C=\chi d,
$$

then

$$
\boxed{
d
=
\Lambda^{1/2}U^\dagger C.
}
$$

The canonical probabilities are

$$
p_k
=
\frac{|d_k|^2}{\sum_l|d_l|^2}.
$$

v0.11 reports the participation ratio

$$
\boxed{
N_{\rm part}
=
\frac{1}{\sum_kp_k^2}.
}
$$

Interpretation:

- $N_{\rm part}\approx1$: one canonical basis direction dominates;
- larger $N_{\rm part}$: the wavefunction genuinely uses multiple independent basis
  directions.

This is more informative than raw $N_{\rm TBF}$ alone.

---

# 24. Spectral effective rank of the overlap matrix

Normalize positive overlap eigenvalues:

$$
w_k
=
\frac{\lambda_k}{\sum_l\lambda_l}.
$$

Define spectral entropy

$$
\boxed{
S_{\rm spec}
=
-\sum_kw_k\ln w_k.
}
$$

The spectral effective rank is

$$
\boxed{
r_{\rm eff}
=
e^{S_{\rm spec}}.
}
$$

This reports how many independent geometric basis directions are represented before
coefficient amplitudes are considered.

---

# 25. v0.11 basis-completeness benchmark

The exact reference observable remains the global-diabatic reduced electronic density

$$
\rho_d^{\rm exact}.
$$

For each Gaussian basis size,

$$
N_{\max}
\in
\{2,4,6,8,10\},
$$

v0.11 computes

$$
\rho_d^{(N_{\max})}.
$$

The principal errors are

$$
\boxed{
\epsilon_P
=
\|
\operatorname{diag}\rho_d^{(N)}
-
\operatorname{diag}\rho_d^{\rm exact}
\|_2,
}
$$

and

$$
\boxed{
\epsilon_\rho
=
\|
\rho_d^{(N)}
-
\rho_d^{\rm exact}
\|_F.
}
$$

The purity error is

$$
\boxed{
\epsilon_{\mathcal P}
=
|
\operatorname{Tr}[(\rho_d^{(N)})^2]
-
\operatorname{Tr}[(\rho_d^{\rm exact})^2]
|.
}
$$

All three are required because matching populations does not guarantee matching
electronic coherence.

---

# 26. Ablation logic

v0.11 runs controlled ablations.

## 26.1 No position optimization

Set

$$
q_c=q_p
$$

for every child.

This tests whether simply changing momentum and width is enough.

## 26.2 Fixed width only

Set

$$
A_c=A_p.
$$

This tests whether width diversity contributes materially.

## 26.3 One child per event

Set

$$
N_{\rm child}=1.
$$

This tests whether multi-child basis growth changes convergence.

The purpose is causal diagnosis.

A feature is considered useful only if removing it measurably worsens the benchmark.

---

# 27. What v0.11 is not

v0.11 is **not** claimed to be a production implementation of the 2009 optimal
spawning algorithm.

Specifically, the present implementation uses:

- a finite local candidate grid rather than continuous nonlinear optimization;
- a first-order local coupling proxy rather than the exact full off-diagonal
  Hamiltonian objective;
- no production AIMS forward/backward spawn-window propagation;
- no ab initio electronic structure in the release benchmark;
- a discrete width bank that is not standard frozen-width AIMS.

The accurate description is:

$$
\boxed{
\text{optimal-spawning-inspired, basis-completeness graph-FMS/AIMS prototype}.
}
$$

That wording is retained throughout the repository.

---

# 28. Scientific success criterion

The v0.11 benchmark is not considered successful merely if the diagonal populations
improve.

The target hierarchy is:

1. generalized norm remains controlled;
2. overlap matrix remains numerically usable;
3. diagonal reduced populations approach exact;
4. reduced-state purity approaches exact;
5. the full reduced density matrix approaches exact;
6. results stabilize under further basis growth.

If populations converge before off-diagonal coherence does, the release reports that
difference explicitly.

That is expected to determine the next scientific bottleneck.
