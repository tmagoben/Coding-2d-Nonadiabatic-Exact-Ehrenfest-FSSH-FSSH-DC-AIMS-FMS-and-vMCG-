# v0.13 Theory: Residual-Driven Gaussian Basis Refinement

Version 0.12 established that, for the strong two-dimensional LVC conical-intersection
benchmark, the dominant remaining error was **not** the projected Gaussian coefficient
propagation.

When the Gaussian method and the exact TDSE started from the same projected initial
state, v0.12 obtained a reduced-density error of only

$$
2.90\times 10^{-4}.
$$

The larger error against the intended exact target came primarily from how well the
coordinate-dependent electronic initial state could be represented by a finite
Gaussian bank.

v0.13 therefore changes the basis-growth question from

> Is the system in a region of strong nonadiabatic coupling?

to

> Which new Gaussian would reduce the unresolved Hilbert-space residual the most?

The release introduces two separate residual concepts:

1. **initial representation residual**, used to construct the initial Gaussian bank;
2. **instantaneous TDSE/Galerkin defect**, used to diagnose missing dynamical tangent
   directions after propagation has begun.

These are related but not interchangeable.

---

# 1. Spinor-complete initial approximation

The v0.12 spinor-complete representation is retained:

$$
\boxed{
\Psi_G(R)
=
\sum_i g_i(R)\mathbf C_i
}
$$

where

$$
\mathbf C_i
=
\begin{pmatrix}
C_{i0}\\
C_{i1}
\end{pmatrix}
$$

is a complete two-state electronic vector in the fixed global diabatic basis.

For a chosen nuclear Gaussian bank

$$
\mathcal B_N
=
\operatorname{span}\{g_1,\ldots,g_N\},
$$

the electronic approximation space is

$$
\boxed{
\mathcal V_N
=
\mathcal B_N\otimes\mathbb{C}^2.
}
$$

The target initial wavefunction is projected orthogonally into this space.

---

# 2. Projection residual

Let

$$
P_N
$$

be the Hilbert-space projector onto

$$
\mathcal V_N.
$$

The current projected state is

$$
\Psi_N=P_N\Psi.
$$

Define the residual

$$
\boxed{
r_N
=
\Psi-\Psi_N.
}
$$

Because this is an orthogonal projection,

$$
\boxed{
\langle v|r_N\rangle=0
\qquad
\forall v\in\mathcal V_N.
}
$$

The squared residual is

$$
\boxed{
\|r_N\|^2
=
\|\Psi- P_N\Psi\|^2.
}
$$

The relative residual reported by the code is

$$
\boxed{
\epsilon_N
=
\frac{\|r_N\|^2}{\|\Psi\|^2}.
}
$$

For normalized targets,

$$
F_N\approx1-\epsilon_N,
$$

where

$$
F_N
$$

is the projection fidelity.

---

# 3. Candidate Gaussian

Consider a normalized new nuclear Gaussian

$$
g_c.
$$

Because every nuclear Gaussian carries the complete two-state electronic basis, the
candidate adds the pair

$$
g_c|d_0\rangle,
\qquad
g_c|d_1\rangle.
$$

However, the part of

$$
g_c
$$

already contained in the existing nuclear span adds no new variational direction.

Define

$$
\boxed{
g_c^\perp
=
(I-P_{\mathcal B_N})g_c.
}
$$

---

# 4. Nonorthogonal formula for the new direction

Let the current nuclear overlap matrix be

$$
S_{ij}=\langle g_i|g_j\rangle.
$$

Define

$$
s_i=\langle g_i|g_c\rangle.
$$

The projection coefficients of the candidate onto the current basis satisfy

$$
S\alpha=s.
$$

Therefore

$$
P_{\mathcal B_N}g_c
=
\sum_i\alpha_i g_i.
$$

Since

$$
\langle g_c|g_c\rangle=1,
$$

the squared norm of the genuinely new direction is

$$
\boxed{
n_c
=
\|g_c^\perp\|^2
=
1-s^\dagger S^{-1}s.
}
$$

If

$$
n_c\rightarrow0,
$$

the candidate is redundant.

This makes residual-driven selection naturally sensitive to basis conditioning.

---

# 5. Exact one-candidate residual gain

Write the electronic residual components as

$$
r_0(R),
\qquad
r_1(R).
$$

Because the current projection residual is orthogonal to the old basis,

$$
\langle g_c^\perp|r_a\rangle
=
\langle g_c|r_a\rangle.
$$

The optimal coefficient added along

$$
g_c^\perp|d_a\rangle
$$

is

$$
\beta_a
=
\frac{
\langle g_c^\perp|r_a\rangle
}{n_c}.
$$

The exact decrease in squared Hilbert residual after adding this Gaussian pair is

$$
\boxed{
\Delta_c
=
\frac{
\displaystyle
\sum_{a=0}^1
|\langle g_c|r_a\rangle|^2
}{n_c}.
}
$$

This is the central v0.13 greedy score.

The pure residual-greedy rule is therefore

$$
\boxed{
c^*
=
\arg\max_c \Delta_c.
}
$$

No nonadiabatic-coupling threshold appears in this criterion.

---

# 6. Why the gain formula is useful

The score answers a very specific question:

> If this Gaussian pair is added right now and the target state is reprojected, how
> much squared wavefunction error can be removed in one step?

This makes the basis expansion directly tied to a measurable approximation error.

For the pure greedy algorithm, the code records both:

- predicted one-step reduction;
- actual residual reduction after the full nonorthogonal reprojection.

They agree within the grid/infinite-domain numerical tolerance tested by the suite.

---

# 7. Deterministic candidate dictionary

v0.13 deliberately avoids a hidden nonlinear optimizer in the release benchmark.

The candidate dictionary is generated from:

```text
center offsets:
    x,y in [-1.0,1.0]
    spacing = 0.2

momentum offsets:
    (0,0)

width scales:
    1.0, 1.5, 2.0, 3.0, 4.0, 6.0
```

This produces

$$
\boxed{726}
$$

deterministic candidate Gaussians.

The dictionary is intentionally transparent.

Changing the dictionary changes the approximation space, so its specification is part
of the numerical method.

---

# 8. Pure residual greedy versus observable-aware screening

The pure greedy algorithm optimizes

$$
\Delta_c.
$$

That is the mathematically cleanest Hilbert-space criterion.

However, v0.12 showed that wavefunction fidelity and a chosen reduced electronic
observable need not improve at exactly the same rate.

Therefore v0.13 includes a second algorithm.

## Stage 1: rigorous residual screen

Rank all admissible candidates by

$$
\Delta_c
$$

and retain only the top

$$
K
$$

candidates.

The release benchmark uses

$$
K=30.
$$

## Stage 2: reduced-density screen

For each of those top residual candidates, perform the one-step trial projection and
compute the initial reduced electronic density error

$$
\epsilon_\rho
=
\|\rho_c-\rho_{\mathrm{target}}\|_F.
$$

Choose the candidate with the smallest

$$
\epsilon_\rho.
$$

Ties are resolved by:

1. smaller wavefunction residual;
2. smaller overlap condition number;
3. larger residual gain.

Thus the observable criterion does **not** replace residual reduction.

It operates only within a set of candidates already certified as strong residual
reducers.

---

# 9. Why density screening is benchmark-specific

The initial exact target state is known in this analytic benchmark.

Therefore its reduced electronic density is available without future information.

The screen uses only

$$
t=0
$$

data.

It does not optimize against the exact final-time answer.

That distinction is essential.

Selecting a basis by inspecting the exact future dynamics would invalidate the
benchmark as a predictive test.

---

# 10. Pure residual result

With 11 Gaussians, pure residual greedy gives

```text
projection fidelity: 0.8910387734744106
relative residual:   0.10896122652558952
initial density error: 0.03719393002837838
condition number:    4664.580844663738
```

The wavefunction residual is lower than in the v0.12 fixed nine-Gaussian bank, but the
initial reduced-density error is not automatically minimized.

That result is scientifically useful:

$$
\boxed{
\text{best Hilbert residual reduction}
\neq
\text{best one-observable reduction at every finite basis size}.
}
$$

---

# 11. Density-screened residual result

The 11-Gaussian v0.13 reference gives

```text
projection fidelity: 0.8902521956060818
relative residual:   0.10974780439391818
initial density error: 0.03209140317550961
condition number:    3458.01502834873
```

The selected sequence remains residual reducing at every step.

The final initial-density error is

$$
\boxed{
\epsilon_{\rho,0}
=
0.032091403
}.
$$

---

# 12. Representation-consistent propagation

Exactly as in v0.12, the final propagation benchmark is not allowed to mix initial
representation error with dynamics error.

The screened v0.13 projected wavefunction is propagated in two ways:

1. exact 2D TDSE;
2. exact-LVC spinor-complete Gaussian propagation.

The projected-state dynamics error is

$$
\boxed{
\epsilon_{\mathrm{dyn}}
=
0.0001135488
}.
$$

This remains much smaller than the target error.

---

# 13. Final target result

Against the intended coordinate-dependent exact benchmark, the v0.13 reference gives

$$
\boxed{
\epsilon_\rho
=
0.031786301
},
$$

population error

$$
\boxed{
\epsilon_P
=
0.025521903
},
$$

trace distance

$$
\boxed{
D_{\mathrm{tr}}
=
0.022476309
},
$$

and coherence phase error

$$
\boxed{
\Delta\phi
=
0.0023799928\;\mathrm{rad}.
}
$$

The generalized norm drift remains

$$
1.059e-06.
$$

---

# 14. v0.12 -> v0.13 comparison

The previous v0.12 fixed nine-Gaussian reference gave

```text
projection fidelity: 0.832276023595292
relative residual: 0.16772397640470793
initial density error: 0.03545457994295867
target density error: 0.03500028070905269
target population error: 0.02810899300694737
condition number: 2235.290713199147
```

v0.13 gives

```text
projection fidelity: 0.8902521956060818
relative residual: 0.10974780439391818
initial density error: 0.03209140317550961
target density error: 0.03178630139393256
target population error: 0.025521902605714804
condition number: 3465.8914579773386
```

The improvement is modest numerically but methodologically important: the v0.13 basis
is generated by a documented error-reduction rule rather than by manually choosing a
fixed grid of Gaussian centers.

---

# 15. Time-dependent Schrödinger defect

Initial projection residual is not enough once propagation begins.

For a moving Gaussian approximation

$$
\Psi_G(t),
$$

define the instantaneous Schrödinger defect

$$
\boxed{
\mathcal R(t)
=
i\frac{d\Psi_G}{dt}
-
\hat H\Psi_G.
}
$$

For an exact solution,

$$
\mathcal R=0.
$$

For a finite Galerkin basis, the projected equations remove the component that can be
represented by the current tangent space, while an unresolved orthogonal component
remains.

v0.13 evaluates this defect explicitly on the exact 2D diagnostic grid.

---

# 16. Coefficient derivative used in the defect

The moving-basis equation is

$$
iS\dot C
=
(H-iT)C.
$$

Therefore

$$
\boxed{
S\dot C
=
-(iH+T)C.
}
$$

At a snapshot, v0.13 solves this linear system for

$$
\dot C.
$$

The total wavefunction derivative is then reconstructed as

$$
\boxed{
\dot\Psi
=
\sum_i
\dot g_i\mathbf C_i
+
\sum_i
g_i\dot{\mathbf C}_i.
}
$$

The same global-diabatic FFT Hamiltonian used by the exact benchmark is applied to
the reconstructed wavefunction.

---

# 17. Galerkin orthogonality diagnostic

For a correctly implemented projected equation, the defect should be nearly
orthogonal to the represented basis directions.

v0.13 therefore reports:

$$
\|\mathcal R\|
$$

and the norm of its projection back into the current spinor-complete basis.

The latter is required to be tiny relative to the total defect in the regression
tests.

This distinguishes a true basis-completeness defect from an error in the projected
coefficient equation.

---

# 18. Residual-driven dynamical enrichment

The same orthogonal-candidate logic can be applied to

$$
\mathcal R.
$$

For a candidate Gaussian pair, define

$$
g_c^\perp
=
(I-P_{\mathcal B})g_c.
$$

The squared TDSE defect that the new pair can capture is

$$
\boxed{
\Delta_c^{\mathrm{TDSE}}
=
\frac{
\displaystyle
\sum_a
|\langle g_c^\perp|\mathcal R_a\rangle|^2
}{\|g_c^\perp\|^2}.
}
$$

The best instantaneous dynamical candidate is

$$
\boxed{
c^*
=
\arg\max_c
\Delta_c^{\mathrm{TDSE}}.
}
$$

---

# 19. Zero-coefficient insertion

The new Gaussian pair enters with

$$
\boxed{
C_{c0}=C_{c1}=0.
}
$$

Therefore the represented wavefunction is exactly unchanged at insertion.

Only the available tangent/Galerkin space changes.

This is a strong audit property:

$$
\boxed{
\Psi_{\mathrm{after\ insertion}}
=
\Psi_{\mathrm{before\ insertion}}.
}
$$

Any immediate change in the physical state would indicate an incorrect basis-growth
implementation.

---

# 20. Actual v0.13 defect-enrichment test

For the release reference, the instantaneous defect enrichment selected

```text
dq=(0.0, -0.4);dp=(0.0, 0.0);width_scale=4
```

The defect norm changed from

$$
0.31502412
$$

to

$$
0.28652498.
$$

The predicted reduction in squared defect was

$$
0.0171436285451,
$$

while the directly recomputed reduction was

$$
0.0171436297853.
$$

The relative prediction error was approximately

$$
7.234e-08.
$$

Thus the residual-capture formula predicts the actual Galerkin defect reduction to
numerical precision for this test.

---

# 21. Why v0.13 does not automatically add this dynamical candidate during the release run

The defect-enrichment primitive is validated, but the primary release propagation
keeps the selected initial bank fixed.

This is deliberate.

A production adaptive algorithm must still answer:

- how often to measure the defect;
- what defect threshold triggers enrichment;
- when a newly added basis function should begin classical guidance;
- whether the child should inherit or optimize a guidance state;
- how to prevent repeated residual additions from degrading conditioning;
- when to remove basis functions.

v0.13 provides the mathematically tested enrichment primitive before embedding it in
a larger state machine.

---

# 22. Relation to AIMS spawning

AIMS/FMS expands a Gaussian basis in nonadiabatic regions, and modern implementations
must control the growth of that basis.

v0.13 is not a replacement for production AIMS spawning.

The residual rule asks a different question:

$$
\boxed{
\text{what unresolved wavefunction direction is missing?}
}
$$

instead of only

$$
\boxed{
\text{where is the nonadiabatic coupling large?}
}
$$

The two criteria can eventually be combined.

For example, a molecular implementation could first identify an electronically
relevant coupling region and then choose among admissible child Gaussians according to
a residual or tangent-space gain.

---

# 23. Relation to variational Gaussian methods

Variational Gaussian methods improve accuracy by allowing the Gaussian basis
parameters themselves to evolve according to variational equations.

v0.13 does not turn the trajectory centers into fully variational degrees of freedom.

Instead it keeps classically guided/frozen TBF propagation and improves the **selection
of basis directions**.

Thus

$$
\boxed{
\text{residual-driven enrichment}
\neq
\text{vMCG parameter evolution}.
}
$$

Both target the finite-basis error, but in different ways.

---

# 24. Conditioning remains part of the algorithm

A candidate with large residual gain can still be numerically dangerous if it is
nearly linearly dependent on the current basis.

v0.13 rejects candidate expansions whose proposed overlap matrix exceeds the configured
condition limit.

The reference selected bank reaches

$$
\kappa(S)_{\max}
\approx
3465.89.
$$

The release acceptance threshold is

$$
5\times10^3.
$$

The separate TDSE-defect enrichment demonstration is allowed to use a looser
diagnostic condition limit because it is not propagated as the release reference.

---

# 25. Scientific conclusion

v0.12 identified the dominant error.

v0.13 turns that diagnosis into an adaptive numerical rule.

The progression is now

$$
\boxed{
\text{coupling-triggered basis growth}
\rightarrow
\text{measured residual-driven basis growth}.
}
$$

The release demonstrates:

1. monotonic Hilbert residual reduction;
2. observable-aware screening without future-time fitting;
3. representation-consistent propagation;
4. explicit TDSE defect evaluation;
5. zero-coefficient defect-driven enrichment;
6. quantitative agreement between predicted and actual defect reduction.

The next logical step is a fully time-adaptive runner that measures the TDSE defect
during propagation and adds/removes TBFs under explicit error and conditioning
budgets.
