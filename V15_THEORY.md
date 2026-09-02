# v0.15 Theory: Shared Gaussian-Pair Algebra and Cost-Aware Adaptation

Version 0.14 established a fully time-adaptive TDSE-defect controller. Its timing audit
identified the dominant implementation cost: repeated unequal-width Gaussian pair
algebra.

v0.15 keeps the v0.14 physical model and changes the numerical architecture. The main
additions are:

1. one reusable Gaussian-pair algebra object per canonical pair;
2. exact conjugate reversal instead of recomputing reversed pair geometry;
3. endpoint pair-cache reuse by the TDSE-defect evaluation;
4. incremental matrix expansion after accepted basis growth;
5. zero-recomputation matrix slicing after pruning;
6. a cost-aware candidate utility that combines predicted defect reduction with the
   incremental computational/stability cost of carrying another TBF.

The release benchmark is deliberately designed so that the physical result remains the
same as v0.14. This isolates the numerical architecture.

## 1. Shared cross-Gaussian algebra

For normalized unequal-width Gaussians, the cross product contains

$$
B_{ij}=A_i+A_j.
$$

Define

$$
\ell_{ij}=A_iq_i+A_jq_j+i(p_j-p_i).
$$

The complex cross centroid and covariance are

$$
\mu_{ij}=B_{ij}^{-1}\ell_{ij},
\qquad
\Sigma_{ij}=B_{ij}^{-1}.
$$

v0.15 obtains both with one multi-right-hand-side solve,

$$
B_{ij}
\begin{pmatrix}
\mu_{ij} & \Sigma_{ij}
\end{pmatrix}
=
\begin{pmatrix}
\ell_{ij} & I
\end{pmatrix}.
$$

The overlap, kinetic energy, LVC potential, and moving-basis time element then reuse
these same quantities.

## 2. Overlap from the cached pair

The overlap is

$$
S_{ij}
=
\exp
\left[
\log P_{ij}
+c_{ij}
+\frac12\ell_{ij}^{T}\mu_{ij}
\right],
$$

with

$$
\log P_{ij}
=
\frac14\log\det A_i
+
\frac14\log\det A_j
-
\frac12\log\det(A_i+A_j)
+
\frac{d}{2}\log 2.
$$

The individual $\log\det A_i$ values are precomputed once per basis snapshot.

## 3. Kinetic and LVC potential from the same pair

With $M^{-1}=B_M$,

$$
u_i=-A_i(\mu_{ij}-q_i)-ip_i,
\qquad
u_j=-A_j(\mu_{ij}-q_j)+ip_j,
$$

and

$$
T_{ij}
=
\frac12 S_{ij}
\left[
u_i^TB_Mu_j
+
\operatorname{Tr}
\left(
A_iB_MA_j\Sigma_{ij}
\right)
\right].
$$

For the two-state LVC potential,

$$
V_d(x,y)
=
\frac12\omega^2(x^2+y^2)I
+
\kappa x\sigma_z
+
\lambda y\sigma_x,
$$

the exact matrix element is

$$
\begin{aligned}
V_{ij}
=
S_{ij}
\Bigg[
&
\frac12\omega^2
\left(
\mu_x^2+\mu_y^2
+\Sigma_{xx}+\Sigma_{yy}
\right)I
\\
&+
\kappa\mu_x\sigma_z
+
\lambda\mu_y\sigma_x
\Bigg].
\end{aligned}
$$

No new pair solve is required.

## 4. Moving-basis time matrix from the same pair

For frozen widths,

$$
T_{ij}^{\mathrm{basis}}
=
\langle g_i|\dot g_j\rangle
=
S_{ij}
\left[
(A_jy-ip_j)^T\dot q_j
+
iy^T\dot p_j
\right],
$$

where

$$
y=\mu_{ij}-q_j.
$$

The general width-derivative terms remain implemented.

## 5. Reverse-pair identity

Because the width matrices are real,

$$
S_{ji}=S_{ij}^*,
\qquad
\mu_{ji}=\mu_{ij}^*,
\qquad
\Sigma_{ji}=\Sigma_{ij}.
$$

Thus one canonical pair $i\le j$ is sufficient even though the moving-basis matrix
itself is not Hermitian.

The cache identity does **not** imply $T_{ji}=T_{ij}^*$. In general,

$$
\dot S=T+T^\dagger.
$$

Only the underlying cross-Gaussian moments are shared.

## 6. v0.14 factorization-equivalent bookkeeping

For one canonical S/H pair, the previous helper structure effectively requested:

```text
outer overlap                              1
kinetic: overlap + centroid + covariance  3
LVC V: overlap + centroid + covariance    3
-------------------------------------------
factorization-equivalent operations       7
```

For one ordered moving-basis pair:

```text
overlap + centroid + covariance = 3
```

These are factorization-equivalent bookkeeping counts, not claims that `solve` and
`inv` have identical low-level constants.

## 7. v0.15 pair snapshot

For one frozen basis snapshot of $N$ Gaussians, v0.15 requires

$$
\frac{N(N+1)}{2}
$$

multi-RHS pair solves.

That one endpoint snapshot can feed:

```text
S
H kinetic blocks
H potential blocks
endpoint TDSE-defect T matrix
```

A separate midpoint snapshot is required because the midpoint Gaussians differ from the
endpoint Gaussians.

At a defect checkpoint, endpoint S/H have already primed every canonical endpoint pair.
Therefore the endpoint defect T matrix requires **zero new pair factorizations**.

## 8. Incremental basis expansion and pruning

If the basis grows $N\rightarrow N+1$ at one fixed endpoint, only the new row/column is
needed. The canonical pair count increases by

$$
\frac{(N+1)(N+2)}{2}
-
\frac{N(N+1)}{2}
=
N+1.
$$

So fixed-snapshot expansion is linear in $N$ pair algebra rather than quadratic.

Moreover, v0.15 constructs a temporary expanded cache while checking the candidate's
condition number. If that candidate is accepted, the same cache is retained. In the
release event, the accepted matrix expansion required

$$
0
$$

new pair factorizations after candidate selection.

Pruning is even cheaper at the same instant: surviving S/H blocks are sliced and the
pair cache is subset/remapped. No pair integral is recomputed.

## 9. Cost-aware adaptation

v0.14 selected the largest TDSE-defect capture fraction. v0.15 first constructs the
same residual-qualified shortlist and then evaluates an incremental computational cost.

The benefit is the predicted fraction of squared defect captured,

$$
f_c
=
\frac{\Delta_c^{\mathrm{TDSE}}}{\|\mathcal R\|^2}.
$$

Over a control horizon of $h$ steps, one extra Gaussian adds $N+1$ canonical pairs to
the endpoint snapshot and $N+1$ to the midpoint snapshot per step:

$$
\Delta P=2h(N+1).
$$

Relative to current pair work, this gives

$$
r_P=\frac{2}{N}.
$$

The dense electronic dimension is $m=sN$, so the relative Cayley overhead is

$$
r_C
=
\frac{[s(N+1)]^3-(sN)^3}{(sN)^3}.
$$

If the horizon contains defect checkpoints, the projected defect solve has the same
cubic ratio $r_D$.

A conditioning penalty is

$$
m_\kappa
=
1+w_\kappa
\max
\left[
\log_{10}
\left(
\frac{\kappa_c}{\kappa_0}
\right),
0
\right].
$$

The release cost is

$$
C_c
=
(r_P+r_C+0.25r_D)m_\kappa,
$$

and utility is

$$
U_c=\frac{f_c}{C_c}.
$$

The release requires

$$
U_c\ge0.15.
$$

The accepted candidate has

$$
U_c=0.24210065,
$$

with normalized incremental cost

$$
C_c=0.62237297.
$$

The ledger also estimates incremental seconds from observed pair-factorization and
Cayley-solve rates. For the release event:

$$
\Delta t_{\mathrm{estimated}}
\approx
0.033310\;\mathrm{s}.
$$

This is a diagnostic extrapolation, not a hard runtime guarantee.

## 10. Release benchmark invariance

The v0.15 benchmark produces

$$
\epsilon_{\mathrm{dyn}}
=
9.5278046e-05,
$$

and

$$
\epsilon_{\mathrm{target}}
=
0.03330494.
$$

The key physical metrics agree with v0.14 to a maximum absolute difference of

$$
8.413e-12.
$$

Thus the cache/incremental architecture changes the implementation, not the benchmark
physics.

## 11. Pair-factorization reduction

The v0.15 propagation performed

```text
actual propagation pair factorizations:
15675

v0.14 factorization-equivalent baseline:
103103
```

for a reduction of

$$
84.80\%.
$$

Candidate conditioning required an additional

```text
88
```

pair factorizations; those are reported separately.

## 12. Runtime comparison

Using the saved v0.14 timing and the v0.15 run in the same build environment:

```text
v0.14 adaptive runtime: 11.289004 s
v0.15 adaptive runtime: 4.207306 s
diagnostic speedup:      2.683 x
runtime reduction:       62.73 %
```

Wall time is environment dependent and is not an acceptance criterion. The portable
claim is the exact removal of repeated pair-factorization work.

## 13. Scientific conclusion

v0.15 is primarily a numerical-architecture release. The physics stays fixed. The gain
comes from recognizing that overlap, kinetic, potential, and moving-basis matrix
elements are different functions of the same cross-Gaussian pair object.

The progression is now

$$
\boxed{
\text{residual-controlled dynamics}
\rightarrow
\text{residual + computational-cost controlled dynamics}.
}
$$

The next scaling step is persistent/local caching across larger basis graphs and sparse
overlap structure once dense $O((sN)^3)$ solves begin to dominate.
