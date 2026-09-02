# v0.10 Theory: Convergence Campaigns, Reduced Electronic Observables, and Error Budgets

Version 0.10 is intentionally different from v0.1-v0.9.

Earlier releases added new dynamical structure.  v0.10 asks whether the accumulated
framework is numerically and scientifically converged.

The central principle is:

$$
\boxed{
\text{a method is not validated by conservation alone;}
\quad
\text{its observables must become insensitive to controlled refinement.}
}
$$

The release therefore separates five questions:

1. Is the exact grid reference converged?
2. Is the Gaussian result converged with timestep?
3. Is it converged with Gaussian basis/spawning controls?
4. Is it insensitive to the SPA0/SPA1 electronic matrix approximation?
5. Does it reproduce the same **well-defined observable** as the exact calculation?

Atomic units are used unless stated otherwise.

---

# 1. Benchmark hierarchy

The benchmark hierarchy is

```text
exact 2D diabatic TDSE
        |
        +-- grid refinement
        +-- timestep refinement
        |
        v
numerical reference
        |
        v
managed graph-Gaussian dynamics
        |
        +-- dt
        +-- SPA order
        +-- spawn action
        +-- overlap blocking
        +-- maximum basis size
        +-- pruning thresholds
        |
        v
observable error + sensitivity budget
```

A Gaussian result is never called "converged" merely because

$$
C^\dagger SC\approx1.
$$

That condition is necessary, but it tests propagation consistency rather than basis or
physical accuracy.

---

# 2. Near-conical-intersection passage benchmark

v0.10 defines a finite-impact-parameter benchmark:

$$
q_0=
\begin{pmatrix}
-0.60\\
0.25
\end{pmatrix},
\qquad
p_0=
\begin{pmatrix}
10\\
0
\end{pmatrix},
$$

with

$$
M=5.
$$

The initial velocity is approximately

$$
\dot q_x(0)=\frac{10}{5}=2.
$$

Over the default interval

$$
t_f=0.60,
$$

the trajectory passes from the negative-$x$ side to the positive-$x$ side of the CI
region while maintaining a nonzero impact parameter in $y$.

The purpose of the nonzero impact parameter is deliberate:

- the trajectory samples strong nonadiabatic coupling;
- the exact adiabatic eigenvectors remain defined along the center path;
- the benchmark avoids placing a numerical TBF center exactly on the singular
  adiabatic point.

The exact wavepacket can still spread around and interfere across the CI region.

---

# 3. Exact-grid convergence

The exact reference uses

$$
e^{-iH\Delta t}
=
e^{-iV\Delta t/2}
e^{-iT\Delta t}
e^{-iV\Delta t/2}
+
\mathcal O(\Delta t^3).
$$

For fixed final time, the global Strang error scales as

$$
\boxed{
\epsilon_t=\mathcal O(\Delta t^2).
}
$$

But timestep convergence does not establish spatial convergence.

The grid spacing is

$$
\Delta x=\frac{2L}{N_x}.
$$

The momentum lattice has spacing

$$
\Delta k=\frac{2\pi}{2L}
$$

and Nyquist magnitude

$$
\boxed{
k_{\max}
\approx
\frac{\pi}{\Delta x}.
}
$$

Therefore a reliable exact reference requires independent control of:

$$
\boxed{
\Delta t,\quad
N_x,\quad
L.
}
$$

v0.10 exposes a two-dimensional `grid_n x dt` surface rather than selecting one
reference calculation without neighboring checks.

---

# 4. Why state-label populations are not sufficient in a graph basis

This is an important correction introduced in v0.10.

Suppose TBF $i$ carries local electronic state label $a_i$ at graph node $n_i$.

It is tempting to define a state population by collecting all coefficients whose
integer label equals $a$.

That is generally **not a projector expectation value** once electronic states live at
different nuclear geometries.

For two TBFs,

$$
\langle
\Phi_{a_i}(R_i)
|
\Phi_{a_j}(R_j)
\rangle
\neq
\delta_{a_i a_j}.
$$

Thus electronic-state sectors defined only by local labels are not mutually
orthogonal.

Consequently, a coefficient-block population can even fail to sum exactly to one.

v0.10 retains the old quantity as a diagnostic proxy for backward compatibility but
does not use it as the primary exact-reference observable.

---

# 5. Reduced electronic density matrix in a common gauge frame

Choose one graph node $r$ as a common electronic reference frame.

Transport each TBF electronic vector to $r$:

$$
\boxed{
|v_i^{(r)}\rangle
=
\mathcal U_{n_i\rightarrow r}
|e_i\rangle.
}
$$

The Gaussian wavefunction in that common frame is approximated as

$$
|\Psi(q)\rangle
=
\sum_i
C_i g_i(q)
|v_i^{(r)}\rangle.
$$

Trace out the nuclear coordinate:

$$
\rho_e^{(r)}
=
\int
|\Psi(q)\rangle
\langle\Psi(q)|
\,dq.
$$

Expanding gives

$$
\boxed{
\rho_e^{(r)}
=
\sum_{ij}
C_i C_j^*
\langle g_j|g_i\rangle
|v_i^{(r)}\rangle
\langle v_j^{(r)}|.
}
$$

Its trace is

$$
\operatorname{Tr}\rho_e^{(r)}
=
C^\dagger SC.
$$

After normalization,

$$
\boxed{
\operatorname{Tr}\rho_e^{(r)}=1.
}
$$

This is a genuine reduced density matrix.

---

# 6. Exact-grid reduced electronic density matrix

The exact diabatic wavefunction is

$$
\Psi_d(R)
=
\begin{pmatrix}
\psi_1(R)\\
\psi_2(R)
\end{pmatrix}.
$$

Tracing over the nuclear grid gives

$$
\boxed{
\rho_d
=
\int
\Psi_d(R)\Psi_d^\dagger(R)
\,dR.
}
$$

Let the common reference-frame electronic matrix be

$$
U_r.
$$

Since

$$
\Psi_d=U_r\Psi_r
$$

for a fixed electronic frame transformation,

$$
\boxed{
\rho_r
=
U_r^\dagger
\rho_d
U_r.
}
$$

Now both the exact and Gaussian reduced electronic density matrices live in the same
fixed electronic basis.

The comparison

$$
\boxed{
\|\rho_{\rm G}-\rho_{\rm exact}\|_F
}
$$

is therefore well defined.

---



# 6A. Preferred analytic-model observable: global diabatic reduced density

The 2D LVC benchmark has an additional advantage that a real molecule usually does
not: it possesses one globally defined diabatic electronic basis.

Therefore the cleanest analytic-model comparison is

$$
\boxed{
\rho_d^{\rm exact}
\quad\text{versus}\quad
\rho_d^{\rm Gaussian}.
}
$$

For Gaussian TBF $i$, the local adiabatic state is converted to its physical vector
in the global diabatic basis,

$$
|u_i^{(d)}\rangle
=
U_{\rm ad}(q_i)|e_{a_i}\rangle.
$$

Then

$$
\boxed{
\rho_d^{\rm Gaussian}
=
\sum_{ij}
C_iC_j^*
S_{ji}^{\rm nuc}
|u_i^{(d)}\rangle
\langle u_j^{(d)}|.
}
$$

The exact density is simply

$$
\boxed{
\rho_d^{\rm exact}
=
\int
\Psi_d(R)
\Psi_d^\dagger(R)
\,dR.
}
$$

This comparison requires no arbitrary reference-node choice.

For the analytic CI benchmark, v0.10 therefore treats the **global diabatic reduced
density** as the preferred exact/Gaussian electronic observable.

For a real ab initio PySCF calculation, where no global diabatic basis is supplied,
the common graph-reference construction remains the appropriate discrete analogue.


# 7. Fixed-frame electronic populations

Once both density matrices are expressed in the same common frame,

$$
\boxed{
P_a^{(r)}
=
(\rho_e^{(r)})_{aa}.
}
$$

These populations satisfy

$$
\sum_aP_a^{(r)}=1.
$$

They are not the same observable as a coordinate-dependent global adiabatic population

$$
\int|\chi_a(R)|^2dR.
$$

v0.10 keeps these two observables explicitly separate.

This prevents a numerically convenient but scientifically invalid comparison between
different definitions of "population."

---

# 8. Electronic purity and electron-nuclear entanglement

For a normalized reduced electronic density matrix,

$$
\boxed{
\mathcal P
=
\operatorname{Tr}(\rho_e^2).
}
$$

For a two-state electronic subsystem,

$$
\frac12\le\mathcal P\le1.
$$

If the full electron-nuclear wavefunction remains separable,

$$
|\Psi\rangle
=
|\chi\rangle|\phi\rangle,
$$

then

$$
\mathcal P=1.
$$

When different nuclear branches become correlated with different electronic states,
tracing out the nuclei produces a mixed electronic state and

$$
\boxed{
\mathcal P<1.
}
$$

The linear entropy is

$$
\boxed{
S_L=1-\mathcal P.
}
$$

The von Neumann entropy is

$$
\boxed{
S_{\rm vN}
=
-\operatorname{Tr}(\rho_e\ln\rho_e)
=
-\sum_k\lambda_k\ln\lambda_k.
}
$$

These diagnostics are especially useful for identifying a Gaussian basis that
preserves norm but fails to generate sufficient nuclear-electronic branching.

---

# 9. Gaussian Wigner initial-condition ensemble

For the frozen Gaussian

$$
g(q)
=
N
\exp
\left[
-\frac12(q-q_0)^TA(q-q_0)
+
ip_0^T(q-q_0)
\right],
$$

with real symmetric positive-definite $A$, the Wigner function is Gaussian in phase
space.

Its coordinate covariance is

$$
\boxed{
\Sigma_q
=
\frac12A^{-1},
}
$$

and its momentum covariance is

$$
\boxed{
\Sigma_p
=
\frac12A.
}
$$

For the zero-chirp packet used here,

$$
\operatorname{Cov}(q,p)=0.
$$

Therefore v0.10 samples

$$
q^{(s)}
\sim
\mathcal N(q_0,\Sigma_q),
$$

$$
p^{(s)}
\sim
\mathcal N(p_0,\Sigma_p)
$$

using a fixed NumPy random seed.

Each sample is then propagated independently.

---

# 10. Ensemble statistics

For observable vector $X_s$ over $N$ sampled initial conditions,

$$
\boxed{
\bar X
=
\frac1N\sum_{s=1}^N X_s.
}
$$

The unbiased sample standard deviation is

$$
\boxed{
\sigma_X
=
\sqrt{
\frac{1}{N-1}
\sum_s
(X_s-\bar X)^2
}.
}
$$

The standard error of the mean is

$$
\boxed{
\operatorname{SEM}
=
\frac{\sigma_X}{\sqrt N}.
}
$$

This does not convert the deterministic quantum problem into a statistical one.

It quantifies sensitivity of the selected initial-condition ensemble and is useful
when the Gaussian method is ultimately used with sampled initial nuclear conditions.

---

# 11. Timestep convergence of the managed Gaussian dynamics

Let

$$
P(\Delta t)
$$

be a final observable.

A three-level refinement uses

$$
\Delta t,
\quad
\frac{\Delta t}{r},
\quad
\frac{\Delta t}{r^2}.
$$

Define successive differences

$$
e_1
=
\|P(\Delta t)-P(\Delta t/r)\|,
$$

$$
e_2
=
\|P(\Delta t/r)-P(\Delta t/r^2)\|.
$$

The observed refinement order is

$$
\boxed{
p_{\rm obs}
=
\frac{\ln(e_1/e_2)}{\ln r}.
}
$$

This estimate is meaningful only when both errors are above roundoff and the
calculation is in an asymptotic refinement regime.

v0.10 therefore records the raw differences as well as the inferred order.

---

# 12. Basis-size convergence

Let

$$
P_{N_G}
$$

be an observable obtained with a maximum of $N_G$ TBFs.

Convergence requires

$$
\boxed{
\|P_{N_G}-P_{N_G+\Delta N}\|
\rightarrow0
}
$$

while maintaining acceptable overlap-matrix conditioning.

A large basis by itself is not evidence of improvement.

If

$$
\kappa(S)\gg1,
$$

the nominally larger basis may contain nearly redundant functions and worsen the
numerical problem.

Thus basis refinement and conditioning must be assessed together.

---

# 13. Repeated spawning

Earlier pedagogical releases prevented a parent TBF from spawning more than once to
the same electronic target.

That makes `max_basis` a poor convergence coordinate for an extended coupling region.

v0.10 adds an **optional repeated-spawning mode**.

After a spawn:

1. the integrated coupling exposure is reset;
2. a minimum number of propagation steps may be required before the same parent-target
   pair is considered again;
3. phase-space overlap blocking still rejects a nearly duplicate child;
4. basis conditioning/pruning remains active.

The default behavior of the inherited v0.9 API remains unchanged unless repeated
spawning is requested.

---

# 14. Overlap-block threshold as a convergence parameter

A proposed child is rejected if its phase-space overlap with an existing target-state
TBF satisfies

$$
\boxed{
|\langle g_{\rm existing}|g_{\rm child}\rangle|
\ge
S_{\rm block}.
}
$$

Changing $S_{\rm block}$ changes how aggressively the basis is allowed to grow.

This is not merely a numerical tolerance.

It is part of the adaptive-basis approximation and must therefore appear in the
convergence surface.

v0.10 includes it as an explicit campaign axis.

---

# 15. SPA-order sensitivity

Let

$$
P^{(0)}
$$

be an observable from SPA0 and

$$
P^{(1)}
$$

the same observable using the first-order electronic Taylor layer.

Define

$$
\boxed{
\epsilon_{\rm SPA}
=
\|P^{(1)}-P^{(0)}\|.
}
$$

A small value indicates local insensitivity to the retained electronic Taylor order.

It does **not** prove all neglected AIMS matrix-element terms are negligible.

That distinction is retained from v0.9.

---

# 16. Spawning-threshold sensitivity

For integrated coupling-action thresholds $A_1$ and $A_2$,

$$
\boxed{
\epsilon_{\rm spawn}
=
\|P(A_1)-P(A_2)\|.
}
$$

A method that changes dramatically under modest threshold refinement is not yet
basis-converged even if each individual run conserves norm.

---

# 17. Multiaxis convergence surface

The managed campaign evaluates the Cartesian product

$$
\boxed{
\{
\Delta t,\;
\text{SPA order},\;
A_{\rm spawn},\;
N_{\rm max},\;
S_{\rm block}
\}.
}
$$

Every row stores:

- final electronic observable;
- generalized norm;
- maximum norm drift;
- maximum overlap condition number;
- maximum basis size;
- spawn count;
- prune count;
- total pruning loss;
- SPA1 matrix correction diagnostic;
- exact-reference error when available.

This makes the convergence study machine-readable rather than a collection of
manually inspected terminal runs.

---

# 18. Exact-reference surface

The exact campaign independently varies

$$
\boxed{
N_{\rm grid}
\quad\text{and}\quad
\Delta t_{\rm exact}.
}
$$

A convenient candidate numerical reference is the row with

$$
\max N_{\rm grid},
\qquad
\min\Delta t.
$$

But v0.10 calls it a **candidate reference** until neighboring refinement results are
checked.

The code therefore names the function `select_finest_exact_reference`, not
`prove_exact_convergence`.

---

# 19. Error and sensitivity budget

v0.10 defines the following population-vector differences:

$$
\epsilon_{\rm exact}
=
\|P_{\rm exact}^{\rm fine}
-
P_{\rm exact}^{\rm next}\|,
$$

$$
\epsilon_{\Delta t}
=
\|P_{\rm G}^{\rm fine}
-
P_{\rm G}^{\rm coarse}\|,
$$

$$
\epsilon_{\rm SPA}
=
\|P_{\rm SPA1}
-
P_{\rm SPA0}\|,
$$

$$
\epsilon_{\rm spawn}
=
\|P_{A_1}
-
P_{A_2}\|,
$$

$$
\epsilon_{\rm basis}
=
\|P_{N_1}
-
P_{N_2}\|,
$$

and the actual reference discrepancy

$$
\boxed{
\epsilon_{\rm total}
=
\|P_{\rm G}-P_{\rm exact}\|.
}
$$

These quantities are **not added in quadrature**.

They are correlated sensitivity probes, not independent random errors.

The largest one identifies the most urgent controlled refinement axis.

---

# 20. Acceptance criteria

A benchmark run can be assigned explicit numerical checks such as

$$
|C^\dagger SC-1|
<
\epsilon_N,
$$

$$
\kappa(S)
<
\kappa_{\max},
$$

$$
\epsilon_{\rm prune}
<
\epsilon_{\rm prune,max},
$$

and, when a validated exact reference is available,

$$
\boxed{
\|P_{\rm G}-P_{\rm exact}\|
<
\epsilon_{\rm ref}.
}
$$

The thresholds are configuration objects in the code.

A failed reference-accuracy check is reported as a failure.

The framework does not reinterpret it as success because conservation tests passed.

---

# 21. PySCF convergence hierarchy

Once the analytic benchmark is understood, the same logic applies to PySCF-driven
direct dynamics.

The electronic layer adds independent convergence axes:

$$
\boxed{
\text{basis set},
\quad
\text{active orbitals},
\quad
\text{active electrons},
\quad
\text{number of states},
\quad
\text{state weights},
\quad
\text{SCF/CASSCF tolerances},
\quad
\text{tracking step}.
}
$$

These must be converged before attributing discrepancies solely to the Gaussian
nuclear method.

v0.10 documents this as a nested hierarchy:

```text
electronic-structure convergence
            ↓
state-tracking/gauge convergence
            ↓
matrix-element approximation convergence
            ↓
Gaussian basis convergence
            ↓
time integration convergence
            ↓
initial-condition ensemble convergence
```

---

# 22. Interpretation philosophy

The most important scientific feature of v0.10 is that a poor comparison is allowed
to remain poor.

For a more demanding passage through the CI region, the exact wavefunction may develop
substantial electronic-nuclear entanglement while a small Gaussian basis remains
nearly pure.

That result is useful.

It identifies which approximation must be improved next.

The purpose of a benchmark framework is not to make every method look accurate.

It is to make the reason for agreement or disagreement measurable.
