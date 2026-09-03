# v0.15 Detailed Derivations

## A. One solve for centroid and covariance

For the cross Gaussian,

$$
B=A_i+A_j,
\qquad
\ell=A_iq_i+A_jq_j+i(p_j-p_i).
$$

Both required quantities satisfy

$$
B\mu=\ell,
\qquad
B\Sigma=I.
$$

With the block right-hand side

$$
R=
\begin{pmatrix}
\ell&I
\end{pmatrix},
$$

one solve gives

$$
BX=R,
\qquad
X=
\begin{pmatrix}
\mu&\Sigma
\end{pmatrix}.
$$

A dense solver can factor $B$ once and apply that factorization to every right-hand
side.

## B. Reverse-pair identities

For the reversed pair,

$$
\ell_{ji}
=
A_jq_j+A_iq_i+i(p_i-p_j)
=
\ell_{ij}^*.
$$

Since $B_{ji}=B_{ij}$ is real,

$$
\mu_{ji}
=
B^{-1}\ell_{ji}
=
\mu_{ij}^*.
$$

Also,

$$
\Sigma_{ji}=\Sigma_{ij},
\qquad
S_{ji}=S_{ij}^*.
$$

So reverse orientation needs no new pair factorization.

## C. Cached kinetic element

With

$$
u_i=-A_i(\mu-q_i)-ip_i,
\qquad
u_j=-A_j(\mu-q_j)+ip_j,
$$

the kinetic element is

$$
K_{ij}
=
\frac12S_{ij}
\left[
u_i^TM^{-1}u_j
+
\operatorname{Tr}
(A_iM^{-1}A_j\Sigma)
\right].
$$

The pair cache supplies $S_{ij}$, $\mu$, and $\Sigma$.

## D. Cached time element

For real time-dependent $A_j$,

$$
\begin{aligned}
\langle g_i|\dot g_j\rangle
=
S_{ij}
\Big[
&
\frac14\operatorname{Tr}(A_j^{-1}\dot A_j)
+
(A_jy-ip_j)^T\dot q_j
+
iy^T\dot p_j
\\
&-
\frac12
\left(
y^T\dot A_jy
+
\operatorname{Tr}(\dot A_j\Sigma)
\right)
\Big],
\end{aligned}
$$

where $y=\mu-q_j$.

For the frozen-width TBFs used by the release, $\dot A_j=0$.

## E. Canonical-pair count

The number of canonical pairs is

$$
P_N
=
\sum_{k=1}^Nk
=
\frac{N(N+1)}{2}.
$$

This is the number of pair-data solves required to prime one snapshot.

## F. v0.14 S/H factorization-equivalent count

For each canonical S/H pair, the old helper structure requested one outer overlap,
three pair operations inside the kinetic helper, and three inside the LVC potential
helper:

$$
F_{SH}^{(0.14)}
=
7P_N
=
7\frac{N(N+1)}{2}.
$$

v0.15 uses

$$
F_{SH}^{(0.15)}
=
P_N.
$$

Therefore the factorization-equivalent reduction for S/H alone is

$$
1-\frac17
=
\frac67
\approx85.714\%.
$$

## G. v0.14 moving-basis count

The old moving-basis helper used overlap, centroid, and covariance for every ordered
pair:

$$
F_T^{(0.14)}=3N^2.
$$

v0.15 uses one canonical midpoint pair snapshot:

$$
F_T^{(0.15)}
=
\frac{N(N+1)}{2}.
$$

Thus

$$
\frac{F_T^{(0.15)}}{F_T^{(0.14)}}
=
\frac{N+1}{6N}.
$$

For large $N$, the reduction approaches $83.33\%$.

## H. Endpoint-defect reuse

At a defect checkpoint, endpoint S/H have already primed the endpoint cache. Hence

$$
F_{T,\mathrm{defect}}^{(0.15)}=0
$$

additional pair factorizations.

The scalar ordered $T_{ij}$ formulas are still evaluated, as are the projected solve,
grid reconstruction, and FFT Hamiltonian.

## I. Incremental addition

Before insertion,

$$
P_N=\frac{N(N+1)}{2}.
$$

After insertion,

$$
P_{N+1}
=
\frac{(N+1)(N+2)}{2}.
$$

The difference is

$$
P_{N+1}-P_N=N+1.
$$

Therefore one fixed-snapshot append requires only $N+1$ new canonical pairs.

## J. Incremental pruning

For retained Gaussian indices $K$, the nuclear overlap is

$$
S_N'=S_N[K,K].
$$

For a spinor-complete basis with $s$ states per Gaussian, let

$$
E=\{sk+a:k\in K,\;a=0,\ldots,s-1\}.
$$

Then

$$
S'=S[E,E],
\qquad
H'=H[E,E].
$$

No integral changes because every surviving basis function is unchanged at the pruning
instant.

## K. Pair cost over a horizon

One extra Gaussian increases both endpoint and midpoint canonical pair snapshots by
$N+1$ each step. Over $h$ steps,

$$
\Delta F_P=2h(N+1).
$$

Current pair work is

$$
2h\frac{N(N+1)}{2}.
$$

Hence

$$
r_P=\frac{2}{N}.
$$

## L. Dense-solve cost

The coefficient dimension is $m=sN$. After adding one TBF,

$$
m'=s(N+1).
$$

The relative cubic work increase is

$$
r_C
=
\frac{[s(N+1)]^3-(sN)^3}{(sN)^3}
=
\frac{(N+1)^3-N^3}{N^3}.
$$

Expanding,

$$
r_C
=
\frac{3}{N}
+
\frac{3}{N^2}
+
\frac{1}{N^3}.
$$

The absolute added cubic work still grows as $O(N^2)$.

## M. Conditioning multiplier

Define

$$
g_\kappa
=
\max
\left[
\log_{10}
\left(
\frac{\kappa_c}{\kappa_0}
\right),
0
\right].
$$

Then

$$
m_\kappa=1+w_\kappa g_\kappa.
$$

This is a numerical-stability penalty, not a physical constant.

## N. Cost-aware utility

With defect capture fraction $f_c$,

$$
C_c
=
(r_P+r_C+\alpha r_D)m_\kappa,
$$

and

$$
U_c=\frac{f_c}{C_c}.
$$

The release uses $\alpha=0.25$.

A candidate must satisfy both a minimum residual capture and a minimum utility. Cost
cannot rescue a candidate that does not reduce the TDSE defect enough.

## O. Empirical seconds estimate

Let $\tau_P$ be observed seconds per pair factorization and $\tau_C$ observed seconds
per dense cubic work unit. The diagnostic incremental-time estimate is

$$
\Delta t_c
=
m_\kappa
\left[
\Delta F_P\tau_P
+
(\Delta W_C+\Delta W_D)\tau_C
\right].
$$

The estimate excludes several smaller grid/control costs and is not a hard timing
guarantee.

## P. Cache memory

A pair object stores a complex overlap, a complex $d$-vector centroid, and a real
$d\times d$ covariance, plus pair-state arrays/references.

The leading mathematical storage is therefore

$$
O(N^2d^2).
$$

For large $d$, dense covariance storage can become significant even before the
coefficient matrix does.

## Q. Physics invariance

The release compares v0.15 with the saved v0.14 reference. The largest difference among
the acceptance metrics is at floating-point level.

The intended invariant is

$$
\boxed{
\text{performance optimization must not alter the represented quantum dynamics}.
}
$$
