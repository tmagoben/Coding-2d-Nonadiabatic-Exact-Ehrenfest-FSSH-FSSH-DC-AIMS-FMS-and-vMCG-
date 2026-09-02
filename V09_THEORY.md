# v0.9 Theory: Convergence-Controlled Graph-AIMS Prototypes

Version 0.9 changes the development goal.

Versions 0.1--0.8 primarily added **new physics and representation machinery**.  v0.9
instead asks whether those algorithms can be made *systematically testable and
convergeable*.

The central hierarchy is

```text
exact 2D TDSE reference
        |
        +--> time-step refinement
        +--> Gaussian-basis refinement
        +--> spawning-threshold refinement
        +--> SPA0 vs first-order Taylor correction
        +--> overlap conditioning / pruning
        |
        v
controlled graph-AIMS benchmark
```

Atomic units are used unless stated otherwise.

---

# 1. Why convergence is a separate scientific problem

A nonadiabatic dynamics program can conserve norm and still be physically inaccurate.
For the present Gaussian framework there are at least five independent errors:

1. **time discretization**;
2. **finite Gaussian basis**;
3. **spawning decision/placement**;
4. **electronic matrix-element approximation**;
5. **electronic-structure/gauge error**.

A useful result must therefore distinguish

$$
\boxed{\text{numerical stability} \neq \text{convergence} \neq \text{chemical accuracy}.}
$$

v0.9 builds explicit diagnostics for the first four while retaining the PySCF/gauge
machinery needed to study the fifth.

---

# 2. Generic Gaussian matrix element

Let two equal-width frozen Gaussians be

$$
g_i(q)
=
N\exp\left[-\frac12(q-q_i)^TA(q-q_i)+ip_i^T(q-q_i)\right],
$$

$$
g_j(q)
=
N\exp\left[-\frac12(q-q_j)^TA(q-q_j)+ip_j^T(q-q_j)\right].
$$

For a smooth scalar or electronic matrix-valued quantity $F(q)$, the central
integral is

$$
\boxed{
M_{ij}[F]
=
\langle g_i|F(q)|g_j\rangle.
}
$$

On-the-fly electronic structure makes direct quadrature of this integral impractical
in many dimensions.  The AIMS saddle-point philosophy replaces the full spatial
dependence by a local Taylor representation near the Gaussian-pair centroid.

---

# 3. Real saddle point for equal widths

Ignoring phase, the magnitude of the Gaussian product is

$$
|g_i(q)g_j(q)|
\propto
\exp\left[
-\frac12(q-q_i)^TA(q-q_i)
-\frac12(q-q_j)^TA(q-q_j)
\right].
$$

Differentiate the exponent:

$$
-A(q-q_i)-A(q-q_j)=0.
$$

Since $A$ is invertible,

$$
2q-q_i-q_j=0.
$$

Therefore the maximum of the product magnitude is

$$
\boxed{
q_c=\frac{q_i+q_j}{2}.
}
$$

This is the real saddle/centroid used by the v0.9 Taylor layer.

---

# 4. Complex cross-Gaussian centroid

The normalized cross density $g_i^*g_j/S_{ij}$ is not a probability density.
Its first moment is complex:

$$
\boxed{
\mu_{ij}
=
\frac{q_i+q_j}{2}
+
\frac{i}{2}A^{-1}(p_j-p_i).
}
$$

Hence

$$
\boxed{
\langle g_i|(q-q_c)|g_j\rangle
=
(\mu_{ij}-q_c)S_{ij}.
}
$$

The imaginary displacement is the reason a first-order Taylor correction need not
vanish even though $q_c$ is the real midpoint.

---

# 5. Zeroth-order saddle-point approximation (SPA0)

Taylor expand a smooth quantity around $q_c$:

$$
F(q)
=
F_c
+
\sum_\alpha F_{c,\alpha}^{(1)}(q_\alpha-q_{c,\alpha})
+
\mathcal O(|q-q_c|^2).
$$

The zeroth-order approximation keeps only

$$
F(q)\approx F_c.
$$

Therefore

$$
\boxed{
M_{ij}^{(0)}[F]
=
F(q_c)S_{ij}.
}
$$

This is the essential structure of the zeroth-order saddle-point approximation used
in practical AIMS: electronic quantities required inside Gaussian integrals are
evaluated at a pair centroid rather than throughout the overlap volume.

---

# 6. First-order Taylor / SPA1 layer

Retaining the linear term gives

$$
M_{ij}^{(1)}[F]
=
F_cS_{ij}
+
\sum_\alpha
F_{c,\alpha}^{(1)}
\langle g_i|(q_\alpha-q_{c,\alpha})|g_j\rangle.
$$

Using the complex first moment,

$$
\boxed{
M_{ij}^{(1)}[F]
=
\left[
F_c
+
\nabla F_c\cdot(\mu_{ij}-q_c)
\right]S_{ij}.
}
$$

For a function exactly linear in the nuclear coordinates, this expression is exact.
The test suite verifies that statement against direct two-dimensional quadrature.

## Terminology caution

The v0.9 code calls this the **SPA1 electronic Taylor layer**.  It is the direct
first-order Taylor extension of the saddle-point electronic matrix element.

It is **not** claimed to reproduce every term of a production AIMS-SPA1 Hamiltonian.
In particular, a complete treatment of first-order expansions of derivative-coupling
terms near conical intersections must be organized so that the final Hamiltonian
remains Hermitian.  The literature discusses this issue explicitly.

---

# 7. Electronic operator form in the graph gauge

At a Gaussian-pair centroid node, v0.8 already supplies a common electronic frame.
Let

$$
H_e(q_c)
$$

be the electronic Hamiltonian and

$$
F_\alpha(q_c)
=
\frac{\partial H_e}{\partial q_\alpha}
$$

its derivative matrices.

Transport the electronic vectors associated with TBFs $i$ and $j$ into that
common frame:

$$
|e_i^{(c)}\rangle,
\qquad
|e_j^{(c)}\rangle.
$$

Define

$$
V_{ij}^{(0)}
=
\langle e_i^{(c)}|H_e(q_c)|e_j^{(c)}\rangle,
$$

and

$$
V_{ij,\alpha}^{(1)}
=
\langle e_i^{(c)}|F_\alpha(q_c)|e_j^{(c)}\rangle.
$$

Then the v0.9 graph-Gaussian potential matrix element is

$$
\boxed{
\langle G_i|V|G_j\rangle_{\mathrm{SPA1}}
=
S_{ij}^{N}
\left[
V_{ij}^{(0)}
+
\sum_\alpha
V_{ij,\alpha}^{(1)}
(\mu_{ij,\alpha}-q_{c,\alpha})
\right].
}
$$

All electronic factors are evaluated in one graph-defined common gauge, so no direct
comparison of unrelated electronic eigenvector phases is required.

---

# 8. Hermiticity of the first-order pair matrix

For equal widths,

$$
S_{ji}=S_{ij}^*,
$$

and

$$
\mu_{ji}=\mu_{ij}^*.
$$

If the electronic operator and derivative matrices are Hermitian and the same
centroid/reference node is used for $(i,j)$ and $(j,i)$, then

$$
V_{ji}^{(0)}=(V_{ij}^{(0)})^*,
$$

$$
V_{ji,\alpha}^{(1)}=(V_{ij,\alpha}^{(1)})^*.
$$

Consequently

$$
\boxed{H_{ji}^{\mathrm{SPA1}}=(H_{ij}^{\mathrm{SPA1}})^*.}
$$

The v0.9 test suite checks this explicitly.

---

# 9. SPA-order difference as an approximation diagnostic

Define

$$
H^{(0)}=H_{\mathrm{SPA0}},
\qquad
H^{(1)}=H_{\mathrm{SPA1}}.
$$

v0.9 reports

$$
\boxed{
\epsilon_{\mathrm{SPA}}
=
\frac{
\|H^{(1)}-H^{(0)}\|_F
}{
\|H^{(1)}\|_F
}.
}
$$

A small $\epsilon_{\mathrm{SPA}}$ means that the first-order electronic variation is
small **for the current Gaussian basis and current overlap regions**.

It does not prove convergence of higher Taylor orders.

---

# 10. Nonorthogonal Gaussian basis conditioning

For

$$
|\Psi\rangle=\sum_iC_i|G_i\rangle,
$$

the overlap matrix is

$$
\boxed{S_{ij}=\langle G_i|G_j\rangle.}
$$

The physical norm is

$$
\boxed{N=C^\dagger SC.}
$$

When two Gaussian basis functions become nearly redundant, $S$ develops a very
small eigenvalue.

Let

$$
S=U\Lambda U^\dagger,
\qquad
\Lambda=\operatorname{diag}(\lambda_1,\ldots,\lambda_N).
$$

The spectral condition number is

$$
\boxed{
\kappa(S)=\frac{\lambda_{\max}}{\lambda_{\min}}.
}
$$

A large $\kappa(S)$ amplifies numerical errors in every solve involving $S$.

---

# 11. Canonical orthogonalization

If eigenvalues below a cutoff are discarded, define

$$
\boxed{
X=U_r\Lambda_r^{-1/2}.
}
$$

Then

$$
\boxed{X^\dagger SX=I.}
$$

v0.9 provides this operation as a diagnostic/regularization utility.

The main dynamic pruning code still retains physical TBFs rather than replacing the
basis by delocalized canonical vectors.

---

# 12. Projecting the wavefunction after deleting a redundant TBF

Suppose the old basis has coefficient vector $C$, and a subset $K$ is retained.
We seek the best projected wavefunction

$$
|\Psi_K\rangle
=
\sum_{i\in K}C_i'|G_i\rangle
$$

that minimizes

$$
\|\Psi-\Psi_K\|^2.
$$

Stationarity gives the normal equations

$$
\boxed{
S_{KK}C'
=
S_{K,\mathrm{all}}C.
}
$$

Thus

$$
\boxed{
C'
=
S_{KK}^{-1}S_{K,\mathrm{all}}C.
}
$$

The old norm is

$$
N_{\mathrm{old}}=C^\dagger SC.
$$

The projected norm is

$$
N_{\mathrm{proj}}=C'^\dagger S_{KK}C'.
$$

Therefore the exact Hilbert-space projection loss is

$$
\boxed{
\epsilon_{\mathrm{prune}}
=
N_{\mathrm{old}}-N_{\mathrm{proj}}\ge0.
}
$$

v0.9 only accepts a deletion if this loss remains below a user-specified budget.

---

# 13. Choosing a redundant TBF

Let $u_{\min}$ be the eigenvector of $S$ associated with the smallest eigenvalue.
It approximately satisfies

$$
\sum_i(u_{\min})_i|G_i\rangle\approx0.
$$

Therefore large components of $u_{\min}$ identify functions participating strongly
in the near-linear dependence.

v0.9 proposes the simplest deterministic policy:

1. find $u_{\min}$;
2. rank basis functions by $|(u_{\min})_i|$;
3. try deleting the largest unprotected candidate;
4. project the wavefunction;
5. accept only if the projection-loss budget is satisfied.

This is intentionally transparent rather than claiming to be the unique optimal AIMS
basis-management algorithm.

---

# 14. Why an instantaneous spawning threshold depends on timestep

The electronic equation contains

$$
\dot c_b
\supset
-
\dot q\cdot d_{ba}\,c_a.
$$

Define

$$
\eta_{ba}(t)
=
|\dot q\cdot d_{ba}|.
$$

Over one sufficiently short time step,

$$
|\Delta c_b|
\sim
\eta_{ba}\Delta t\,|c_a|.
$$

Therefore a criterion formulated solely as

$$
\eta>\eta_{\mathrm{threshold}}
$$

is not directly tied to the amplitude transferred during the numerical step.

---

# 15. Integrated coupling-action spawning criterion

v0.9 accumulates

$$
\boxed{
\mathcal A_{ab}(t)
=
\int_{\mathrm{coupling\ region}}
|\dot q\cdot d_{ab}|\,dt.
}
$$

Numerically,

$$
\boxed{
\mathcal A_{ab}^{n+1}
=
\mathcal A_{ab}^{n}
+
|\dot q\cdot d_{ab}|_n\Delta t.
}
$$

A child becomes eligible when

$$
\boxed{
\mathcal A_{ab}>\mathcal A_{\mathrm{spawn}}.
}
$$

If the instantaneous coupling falls below a small floor, the coupling-region
accumulator is reset.

This is a v0.9 **adaptive prototype criterion**, not a claim that conventional AIMS
uses this exact integrated threshold.

Its advantage is clear: for a constant coupling rate, the trigger time converges as
$\Delta t$ is refined.

---

# 16. Child placement remains energy constrained

The v0.8 NAC-direction child momentum is retained.

For mass metric $B=M^{-1}$,

$$
p_b=p_a+\lambda n,
$$

with $n$ along the derivative-coupling direction.

Energy conservation requires

$$
\boxed{
(n^TBn)\lambda^2
+2(p_a^TBn)\lambda
+2(E_b-E_a)=0.
}
$$

A negative discriminant means the current local placement rule cannot construct a
real energy-conserving child.

---

# 17. Time-step convergence

Let a numerical observable be $Q(\Delta t)$.
For a method of global order $p$,

$$
Q(h)=Q^*+Ch^p+\mathcal O(h^{p+1}).
$$

For refinements $h,h/r,h/r^2$, define successive differences

$$
e_h=\|Q(h)-Q(h/r)\|,
$$

$$
e_{h/r}=\|Q(h/r)-Q(h/r^2)\|.
$$

Then the observed order is

$$
\boxed{
p_{\mathrm{obs}}
=
\frac{\ln(e_h/e_{h/r})}{\ln r}.
}
$$

A large or erratic $p_{\mathrm{obs}}$ can occur when the observable is already at
roundoff-level convergence; observed order should therefore be interpreted together
with the absolute error scale.

---

# 18. Basis-size convergence

A spawned Gaussian calculation has another refinement axis:

$$
N_{\mathrm{TBF}}=1,2,3,\ldots.
$$

The relevant comparison is not just coefficient-vector dimension because the basis
changes.

Use physical observables such as

$$
P_I(t),
\qquad
\langle q\rangle,
\qquad
\text{final channel yield},
$$

and monitor

$$
\boxed{
\epsilon_N
=
\|Q_{N_{\mathrm{TBF}}}-Q_{N_{\mathrm{TBF}}+\Delta N}\|.
}
$$

A calculation whose answer changes substantially when one additional allowed spawn is
included is not basis converged.

---

# 19. Spawning-threshold convergence

Likewise vary

$$
\mathcal A_{\mathrm{spawn}}.
$$

A practical convergence sequence is

$$
\mathcal A,
\quad
\mathcal A/2,
\quad
\mathcal A/4.
$$

Record

- number of TBFs;
- spawn times;
- final populations;
- maximum $\kappa(S)$;
- pruning events;
- electronic-structure evaluations.

A threshold is not converged merely because two runs happened to spawn the same
number of children.

---

# 20. Exact two-dimensional TDSE reference

The reference Hamiltonian remains

$$
\boxed{
H_d
=
-\frac{1}{2M}(\partial_x^2+\partial_y^2)I_2
+V_d(x,y).
}
$$

It is propagated with second-order Strang splitting:

$$
\boxed{
e^{-iH\Delta t}
\approx
e^{-iV\Delta t/2}
e^{-iT\Delta t}
e^{-iV\Delta t/2}.}
$$

The exact-grid calculation is still numerical; it must itself be converged with
respect to

$$
\Delta x,\Delta y,\Delta t,
$$

and box size.

Within a sufficiently converged grid, however, it provides the most useful available
reference for the analytic two-dimensional CI model.

---

# 21. Matching the exact and Gaussian initial states

v0.9 initializes the exact wavefunction as

$$
\boxed{
\Psi_d(R,0)
=
g(R;q_0,p_0,A)\Phi_a(R),
}
$$

where $\Phi_a(R)$ is the selected adiabatic electronic state represented in the
diabatic basis.

This mirrors the interpretation of one initial adiabatic Gaussian TBF.

Because a global real adiabatic gauge around a CI is impossible, the analytic frame
contains the expected branch/sign structure.  Population observables remain gauge
invariant.

---

# 22. Exact adiabatic populations

At every grid point,

$$
\chi(R)=U^\dagger(R)\Psi_d(R).
$$

Then

$$
\boxed{
P_a
=
\int|\chi_a(R)|^2dR.
}
$$

The test suite verifies

$$
\boxed{P_0+P_1=\|\Psi\|^2}
$$

to numerical precision.

---

# 23. Managed graph-AIMS workflow

The new `run_managed_graph_aims` propagator combines the preceding pieces:

```text
1. build graph + TBF Gaussian matrices
2. choose SPA0 or SPA1 electronic Taylor order
3. move TBF centers
4. update temporal and centroid graph nodes
5. construct metric-compatible moving-basis T
6. propagate i S Cdot = (H - iT) C
7. accumulate nonadiabatic coupling action
8. spawn when the action budget is exceeded
9. rebuild S,H after spawning
10. monitor cond(S)
11. if necessary, projectively prune a redundant TBF
12. record norm, populations, graph size, SPA correction, events
```

Every approximation is therefore represented by an explicit setting and diagnostic.

---

# 24. Three levels of validation

## Level I: algebraic invariants

Examples:

$$
S=S^\dagger,
\qquad
H=H^\dagger,
\qquad
C^\dagger SC\approx1.
$$

These detect coding errors.

## Level II: numerical self-convergence

Examples:

$$
\Delta t\rightarrow\Delta t/2,
\qquad
N_{\mathrm{TBF}}\rightarrow N_{\mathrm{TBF}}+1.
$$

These detect discretization and truncation error.

## Level III: external/reference convergence

Compare against the exact 2D TDSE:

$$
\boxed{
\epsilon_P
=
\|P^{\mathrm{Gaussian}}-P^{\mathrm{exact}}\|_2.
}
$$

This tests the physical approximation itself.

Passing Level I is necessary but never sufficient.

---

# 25. What v0.9 calls "converged"

A v0.9 benchmark should only be labeled numerically converged after reporting, at
minimum:

1. exact-grid refinement;
2. Gaussian timestep refinement;
3. spawning-action refinement;
4. maximum-basis refinement;
5. SPA0/SPA1 sensitivity;
6. maximum $\kappa(S)$;
7. cumulative pruning projection loss;
8. norm drift;
9. final population error against the exact reference.

The repository provides the machinery for these checks but does not hard-code one
universal tolerance for every molecular problem.

---

# 26. Relation to conventional AIMS

The AIMS literature uses the saddle-point approximation to replace expensive nuclear
integrals of electronic quantities by evaluations near TBF-pair centroids.  Practical
AIMS commonly employs a zeroth-order saddle-point approximation, together with other
approximations such as the independent first generation approximation.

v0.9 deliberately separates three ideas:

- **SPA0-like centroid evaluation**;
- **first-order electronic Taylor correction**;
- **graph/overlap gauge transport**.

This makes each approximation independently testable on the analytic CI model.

It should not be presented as a bit-for-bit reimplementation of a production AIMS
package.

---

# 27. Remaining limitations

Even after v0.9, a production molecular AIMS/vMCG program would still need more mature
implementations of:

- complete AIMS SPA0/SPA1 interstate/intrastate matrix elements;
- optimized spawning;
- IFGA and alternatives handled explicitly;
- multidimensional ab initio initial-condition ensembles;
- adaptive electronic-structure failure recovery;
- large-scale TBF group splitting/pruning;
- nonlinear molecular coordinates and mass metrics;
- statistical convergence over initial conditions;
- extensive molecule-specific active-space and state-manifold convergence.

v0.9's goal is more fundamental:

$$
\boxed{
\text{make approximation error visible before increasing chemical complexity.}
}
$$
