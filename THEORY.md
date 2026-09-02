# Theory: Gaussian Wavepacket Dynamics from Heller to Multiple Spawning

This document develops the framework in the same order as the code.

Atomic units are used throughout unless explicitly stated:

$$
\hbar = 1.
$$

The one-dimensional Hamiltonian is

$$
\boxed{
\hat H
=
-\frac{1}{2m}\frac{\partial^2}{\partial x^2}
+
V(x).
}
$$

The time-dependent Schrödinger equation (TDSE) is

$$
\boxed{
i\frac{\partial}{\partial t}\Psi(x,t)
=
\hat H\Psi(x,t).
}
$$

The central question is:

> How much of the exact wavefunction can be represented by one or a small number of
> moving Gaussian functions?

---

# 1. A normalized Gaussian wavepacket

We begin with

$$
g(x;q,p,\alpha)
=
N_\alpha
\exp\left[
-\frac{\alpha}{2}(x-q)^2
+
ip(x-q)
\right],
$$

where

$$
\alpha>0.
$$

The parameters have direct meanings:

- $q$: position center,
- $p$: mean momentum,
- $\alpha$: inverse-square width.

The normalization constant follows from

$$
1
=
\int_{-\infty}^{\infty}|g(x)|^2 dx
=
|N_\alpha|^2
\int_{-\infty}^{\infty}
e^{-\alpha(x-q)^2}dx.
$$

Using

$$
\int_{-\infty}^{\infty}
e^{-\alpha y^2}dy
=
\sqrt{\frac{\pi}{\alpha}},
$$

we obtain

$$
\boxed{
N_\alpha
=
\left(\frac{\alpha}{\pi}\right)^{1/4}.
}
$$

The coordinate probability density is therefore

$$
|g(x)|^2
=
\sqrt{\frac{\alpha}{\pi}}
e^{-\alpha(x-q)^2}.
$$

Comparing this with the normal-distribution form gives

$$
\boxed{
\langle x\rangle=q,
\qquad
\operatorname{Var}(x)=\frac{1}{2\alpha}.
}
$$

Thus the standard deviation is

$$
\sigma_x=\frac{1}{\sqrt{2\alpha}}.
$$

---

# 2. Momentum expectation and uncertainty

Differentiate the Gaussian:

$$
\frac{\partial g}{\partial x}
=
[-\alpha(x-q)+ip]g.
$$

Therefore,

$$
\hat p g
=
-i\frac{\partial g}{\partial x}
=
[p+i\alpha(x-q)]g.
$$

The odd term integrates to zero, so

$$
\boxed{
\langle \hat p\rangle=p.
}
$$

A second derivative gives

$$
\frac{\partial^2 g}{\partial x^2}
=
\left(
[-\alpha(x-q)+ip]^2-\alpha
\right)g.
$$

From this one obtains

$$
\operatorname{Var}(p)=\frac{\alpha}{2}.
$$

Hence

$$
\boxed{
\sigma_x\sigma_p
=
\frac12.
}
$$

The real-width Gaussian is therefore a minimum-uncertainty state.

---

# 3. Heller's complex Gaussian ansatz

To permit spreading and position-momentum correlation, replace the fixed real width
by a complex width parameter.

Use

$$
\boxed{
\Psi(x,t)
=
\exp
\left\{
i\left[
\frac12 A_t(x-q_t)^2
+
p_t(x-q_t)
+
\gamma_t
\right]
\right\}.
}
$$

Here

$$
A_t=A_R(t)+iA_I(t),
$$

with

$$
A_I(t)>0
$$

for a square-integrable packet.

The Gaussian envelope is determined by

$$
\left|
e^{iA(x-q)^2/2}
\right|
=
e^{-A_I(x-q)^2/2}.
$$

Thus

$$
\boxed{
\sigma_x^2=\frac{1}{2A_I}.
}
$$

The real part $A_R$ is a quadratic phase, or "chirp", which represents
position-momentum correlation.

The scalar $\gamma_t$ is complex. Its real part contains dynamical phase information
and its imaginary part controls normalization.

---

# 4. Local harmonic expansion of the potential

Heller's thawed Gaussian approximation begins by expanding the potential around the
instantaneous packet center:

$$
V(x)
=
V(q)
+
V'(q)(x-q)
+
\frac12V''(q)(x-q)^2
+
\mathcal O[(x-q)^3].
$$

The TGA discards cubic and higher terms:

$$
\boxed{
V(x)
\approx
V(q)
+
V'(q)\xi
+
\frac12V''(q)\xi^2,
\qquad
\xi=x-q.
}
$$

This is the essential approximation.

If the physical potential is globally quadratic, the expansion is exact and the
single Gaussian remains an exact Gaussian solution of the TDSE.

---

# 5. Derivation of the thawed-Gaussian equations

Define the complex phase polynomial

$$
S(x,t)
=
\frac12 A\xi^2+p\xi+\gamma,
\qquad
\xi=x-q.
$$

Then

$$
\Psi=e^{iS}.
$$

## 5.1 Time derivative

Because

$$
\dot\xi=-\dot q,
$$

we have

$$
\dot S
=
\frac12\dot A\xi^2
-
A\xi\dot q
+
\dot p\,\xi
-
p\dot q
+
\dot\gamma.
$$

Since

$$
i\partial_t\Psi
=
-\dot S\,\Psi,
$$

the left side of the TDSE is

$$
-\dot S.
$$

## 5.2 Spatial derivatives

The first spatial derivative is

$$
\frac{\partial S}{\partial x}
=
A\xi+p.
$$

The second derivative is

$$
\frac{\partial^2 S}{\partial x^2}=A.
$$

For $\Psi=e^{iS}$,

$$
-\frac{1}{2m}
\frac{\partial^2\Psi}{\partial x^2}
=
\left[
\frac{(A\xi+p)^2}{2m}
-
\frac{iA}{2m}
\right]\Psi.
$$

## 5.3 Match powers of $\xi$

Substitute the local quadratic potential and collect equal powers of $\xi$.

### Quadratic term

$$
-\frac12\dot A
=
\frac{A^2}{2m}
+
\frac12V''(q),
$$

so

$$
\boxed{
\dot A
=
-\frac{A^2}{m}
-
V''(q).
}
$$

This Riccati equation controls spreading and chirp.

### Linear term

The left-hand coefficient is

$$
A\dot q-\dot p.
$$

The right-hand coefficient is

$$
\frac{Ap}{m}+V'(q).
$$

Choosing

$$
\boxed{
\dot q=\frac{p}{m}
}
$$

leaves

$$
\boxed{
\dot p=-V'(q).
}
$$

Thus the Gaussian center follows classical Hamiltonian mechanics.

### Constant term

Using $\dot q=p/m$,

$$
\frac{p^2}{m}-\dot\gamma
=
\frac{p^2}{2m}
+
V(q)
-
\frac{iA}{2m}.
$$

Hence

$$
\boxed{
\dot\gamma
=
\frac{p^2}{2m}
-
V(q)
+
\frac{iA}{2m}.
}
$$

Define the classical Lagrangian

$$
L
=
\frac{p^2}{2m}-V(q).
$$

Then

$$
\boxed{
\dot\gamma
=
L+\frac{iA}{2m}.
}
$$

The real dynamical phase therefore contains the classical action, while the complex
width contributes the quantum normalization/focusing term.

---

# 6. Heller TGA equations collected

The complete one-dimensional TGA system is

$$
\boxed{
\dot q=\frac{p}{m},
}
$$

$$
\boxed{
\dot p=-V'(q),
}
$$

$$
\boxed{
\dot A=-\frac{A^2}{m}-V''(q),
}
$$

$$
\boxed{
\dot\gamma
=
\frac{p^2}{2m}
-
V(q)
+
\frac{iA}{2m}.
}
$$

These are exactly the equations implemented in `gaussian_dynamics/heller.py`.

---

# 7. Initializing the complex width

Suppose the desired initial normalized packet is

$$
\Psi(x,0)
=
\left(\frac{1}{2\pi\sigma^2}\right)^{1/4}
\exp\left[
-\frac{(x-q_0)^2}{4\sigma^2}
+
ip_0(x-q_0)
\right].
$$

Compare the quadratic terms:

$$
\frac{i}{2}A_0(x-q_0)^2
=
-\frac{(x-q_0)^2}{4\sigma^2}.
$$

Therefore,

$$
\boxed{
A_0=\frac{i}{2\sigma^2}.
}
$$

The initial normalization factor is

$$
N_0
=
(2\pi\sigma^2)^{-1/4}.
$$

Since

$$
e^{i\gamma_0}=N_0,
$$

we may choose

$$
\boxed{
\gamma_0=-i\ln N_0.
}
$$

---

# 8. Why TGA is exact for quadratic potentials

Let

$$
V(x)=V_0+f x+\frac12 kx^2.
$$

A Taylor expansion around any point $q$ terminates exactly at second order because

$$
V'''(x)=0.
$$

Therefore no discarded term exists.

The TDSE acting on a Gaussian produces at most a quadratic polynomial times that
Gaussian, which lies in the tangent space generated by changes in

$$
q,\;p,\;A,\;\gamma.
$$

Thus the Gaussian manifold is closed under quadratic Hamiltonian evolution.

This yields the important result

$$
\boxed{
\text{Heller TGA is exact for a single Gaussian under any quadratic Hamiltonian.}
}
$$

This is explicitly tested against the split-operator reference.

---

# 9. Frozen Gaussian approximation

A simpler trajectory-guided ansatz fixes the real width:

$$
\boxed{
g_F(x,t)
=
\left(\frac{\alpha}{\pi}\right)^{1/4}
\exp
\left[
-\frac{\alpha}{2}(x-q_t)^2
+
ip_t(x-q_t)
+
iS_t
\right].
}
$$

The center follows

$$
\dot q=\frac{p}{m},
\qquad
\dot p=-V'(q),
$$

and the trajectory action obeys

$$
\dot S=\frac{p^2}{2m}-V(q).
$$

This ansatz preserves its width by construction.

It is less flexible than the TGA because it cannot represent breathing/spreading with
one basis function.

Its importance is conceptual: **frozen Gaussians are the building blocks used in
multiple-spawning approaches such as FMS/AIMS.**

---

# 10. Why one Gaussian eventually fails

For a strongly anharmonic potential, the exact wavepacket can:

- distort,
- skew,
- split,
- develop nodes,
- interfere with itself,
- tunnel.

A single Gaussian has only a finite number of shape parameters. No choice of
$q,p,A,\gamma$ can represent a two-peaked wavefunction exactly.

The natural next step is therefore

$$
\boxed{
\Psi(x,t)
=
\sum_{j=1}^{N_G}
C_j(t)g_j(x,t).
}
$$

---

# 11. A nonorthogonal moving Gaussian basis

Let

$$
|\Psi\rangle
=
\sum_j C_j(t)|g_j(t)\rangle.
$$

Differentiate:

$$
|\dot\Psi\rangle
=
\sum_j \dot C_j|g_j\rangle
+
\sum_j C_j|\dot g_j\rangle.
$$

Insert into the TDSE,

$$
i|\dot\Psi\rangle=\hat H|\Psi\rangle.
$$

Project with $\langle g_i|$:

$$
i
\sum_j
\langle g_i|g_j\rangle\dot C_j
+
i
\sum_j
\langle g_i|\dot g_j\rangle C_j
=
\sum_j
\langle g_i|\hat H|g_j\rangle C_j.
$$

Define

$$
\boxed{
S_{ij}=\langle g_i|g_j\rangle,
}
$$

$$
\boxed{
H_{ij}=\langle g_i|\hat H|g_j\rangle,
}
$$

$$
\boxed{
\tau_{ij}=\langle g_i|\dot g_j\rangle.
}
$$

Then

$$
\boxed{
iS\dot C=(H-i\tau)C.
}
$$

Equivalently,

$$
\boxed{
\dot C
=
-iS^{-1}(H-i\tau)C.
}
$$

For nearly linearly dependent Gaussians, explicitly forming $S^{-1}$ is numerically
poor practice. The code solves the linear system directly.

---

# 12. Why the $\tau$ matrix cannot be omitted

If the basis depends on time, changing the Gaussian basis itself contributes to
$\dot\Psi$.

Omitting $\tau$ incorrectly treats the moving basis as stationary.

Differentiate the overlap:

$$
\dot S_{ij}
=
\langle\dot g_i|g_j\rangle
+
\langle g_i|\dot g_j\rangle.
$$

Therefore,

$$
\boxed{
\dot S=\tau^\dagger+\tau.
}
$$

This identity is tested numerically.

It is also the key relation that permits norm conservation in an exactly projected
time-dependent basis.

The wavefunction norm is

$$
\boxed{
\langle\Psi|\Psi\rangle=C^\dagger S C.
}
$$

One must monitor this quantity, not merely $C^\dagger C$, because the basis is
nonorthogonal.

---

# 13. Frozen Gaussian basis derivatives

For

$$
g_j(x)
=
N
\exp
\left[
-\frac{\alpha}{2}(x-q_j)^2
+
ip_j(x-q_j)
\right],
$$

with fixed $\alpha$,

$$
\frac{\partial g_j}{\partial q_j}
=
\left[
\alpha(x-q_j)-ip_j
\right]g_j,
$$

and

$$
\frac{\partial g_j}{\partial p_j}
=
i(x-q_j)g_j.
$$

Therefore,

$$
\boxed{
\dot g_j
=
\frac{\partial g_j}{\partial q_j}\dot q_j
+
\frac{\partial g_j}{\partial p_j}\dot p_j.
}
$$

If the centers are trajectory guided,

$$
\dot q_j=\frac{p_j}{m},
\qquad
\dot p_j=-V'(q_j).
$$

These derivatives are used directly in the moving-basis implementation.

---

# 14. Analytic overlap of equal-width frozen Gaussians

For equal real width $\alpha$,

$$
g_j(x)
=
\left(\frac{\alpha}{\pi}\right)^{1/4}
e^{-\alpha(x-q_j)^2/2+ip_j(x-q_j)},
$$

the overlap can be integrated analytically:

$$
\boxed{
S_{ij}
=
\exp
\left[
-\frac{\alpha}{4}(q_i-q_j)^2
-
\frac{(p_i-p_j)^2}{4\alpha}
+
\frac{i}{2}(p_i+p_j)(q_i-q_j)
\right].
}
$$

This formula has an immediate phase-space interpretation.

Two Gaussians become nearly orthogonal when they are separated strongly in position,

$$
|q_i-q_j|\gg \alpha^{-1/2},
$$

or momentum,

$$
|p_i-p_j|\gg \alpha^{1/2}.
$$

The repository compares this analytic expression with grid quadrature.

---

# 15. Moving-basis MCG foundation

A simple multiple-Gaussian construction uses:

1. classical equations for each Gaussian center;
2. fixed Gaussian widths;
3. quantum coefficients propagated through

$$
iS\dot C=(H-i\tau)C.
$$

This already differs fundamentally from an ensemble of independent trajectories.

The full wavefunction is coherent:

$$
|\Psi|^2
=
\sum_{ij}
C_i^*C_j
g_i^*g_j.
$$

The terms with $i\ne j$ are interference terms.

Thus two Gaussian basis functions can interfere even if their centers follow
classical-looking paths.

This repository calls this layer the **moving-basis MCG foundation** because the
Gaussian centers are prescribed rather than fully determined by a variational
principle.

---

# 16. Time-dependent variational principle

The next step is to allow the Gaussian parameters themselves to respond to the
quantum wavefunction.

Let

$$
|\Psi(\boldsymbol\theta)\rangle
$$

be a variational wavefunction depending on real parameters

$$
\theta_1,\ldots,\theta_M.
$$

Then

$$
|\dot\Psi\rangle
=
\sum_\nu
\left|
\frac{\partial\Psi}{\partial\theta_\nu}
\right\rangle
\dot\theta_\nu.
$$

Define tangent vectors

$$
|D_\nu\rangle
=
\left|
\frac{\partial\Psi}{\partial\theta_\nu}
\right\rangle.
$$

---

# 17. McLachlan variational principle

Define the TDSE residual

$$
|r\rangle
=
i|\dot\Psi\rangle-\hat H|\Psi\rangle.
$$

McLachlan's principle chooses $\dot{\boldsymbol\theta}$ to minimize

$$
\boxed{
\|r\|^2.
}
$$

Substitute

$$
|\dot\Psi\rangle
=
\sum_\nu |D_\nu\rangle\dot\theta_\nu.
$$

Then

$$
|r\rangle
=
i\sum_\nu |D_\nu\rangle\dot\theta_\nu
-
\hat H|\Psi\rangle.
$$

Differentiate $\langle r|r\rangle$ with respect to the real velocity
$\dot\theta_\mu$ and set the result to zero:

$$
\operatorname{Re}
\langle iD_\mu|r\rangle
=
0.
$$

This gives

$$
\sum_\nu
\operatorname{Re}
\langle D_\mu|D_\nu\rangle
\dot\theta_\nu
=
\operatorname{Im}
\langle D_\mu|\hat H|\Psi\rangle.
$$

Define

$$
\boxed{
G_{\mu\nu}
=
\operatorname{Re}
\langle D_\mu|D_\nu\rangle,
}
$$

and

$$
\boxed{
b_\mu
=
\operatorname{Im}
\langle D_\mu|\hat H|\Psi\rangle.
}
$$

The variational equations are

$$
\boxed{
G\dot{\boldsymbol\theta}=b.
}
$$

This equation is the core of `gaussian_dynamics/variational.py`.

---

# 18. Why this is useful for understanding vMCG

In a production vMCG derivation, one usually adopts carefully chosen complex Gaussian
parameters and derives coupled coefficient/parameter equations analytically, often
within a G-MCTDH formalism.

For a pedagogical implementation, the tangent-space equation above is valuable because
it makes the variational idea explicit:

1. choose a multi-Gaussian manifold;
2. calculate the tangent vectors;
3. project the exact TDSE velocity into that tangent space;
4. solve a small linear system;
5. propagate the variational parameters.

The code uses numerical derivatives of the wavefunction with respect to the parameters.
This is slower than analytic derivatives but substantially easier to audit line by
line.

The framework therefore teaches the **variational mechanism** without pretending to
be a production DD-vMCG implementation.

---

# 19. Parameterization used in the variational example

For each Gaussian $j$, use

$$
g_j(x)
=
\left(\frac{\alpha_j}{\pi}\right)^{1/4}
\exp
\left[
-\frac{\alpha_j}{2}(x-q_j)^2
+
ip_j(x-q_j)
\right].
$$

The total wavefunction is

$$
\Psi(x)
=
\sum_j
(c_{j,R}+ic_{j,I})g_j(x).
$$

To ensure

$$
\alpha_j>0,
$$

write

$$
\boxed{
\alpha_j=e^{\lambda_j}.
}
$$

To permit width-position/momentum correlation, include a real quadratic phase
(chirp) $\kappa_j$:

$$
g_j(x)
=
\left(\frac{\alpha_j}{\pi}\right)^{1/4}
\exp
\left[
-\frac{\alpha_j}{2}(x-q_j)^2
+
ip_j(x-q_j)
+
\frac{i\kappa_j}{2}(x-q_j)^2
\right].
$$

The real variational parameter set is therefore

$$
\boxed{
\theta
=
\{
c_{j,R},c_{j,I},q_j,p_j,\lambda_j,\kappa_j
\}_{j=1}^{N_G}.
}
$$

The pair $(\alpha_j,\kappa_j)$ is the real-parameter counterpart of a complex
Gaussian width. This is not the most efficient parameterization, but every parameter
is transparent.

---

# 20. Linear dependence and regularization

If two Gaussians become nearly identical, their tangent vectors become nearly linearly
dependent.

Then $G$ can become ill conditioned.

The code therefore uses a least-squares/pseudoinverse solve with an explicit
singular-value cutoff rather than blindly evaluating

$$
G^{-1}.
$$

This is not merely a software issue. It reflects redundancy of the variational
coordinates.

A serious high-dimensional vMCG implementation requires more sophisticated handling
of these conditioning problems.

---

# 21. Multielectronic-state Born-Huang expansion

To connect Gaussian dynamics with nonadiabatic dynamics, return to the molecular
wavefunction:

$$
\boxed{
\Psi(r,R,t)
=
\sum_I
\Omega_I(R,t)
\Phi_I(r;R).
}
$$

Here

- $I$ labels electronic states,
- $\Phi_I$ are adiabatic electronic states,
- $\Omega_I$ are nuclear wavefunctions.

A Gaussian-basis representation expands each nuclear wavefunction:

$$
\boxed{
\Omega_I(R,t)
=
\sum_{k=1}^{N_I}
C_k^{(I)}(t)
\chi_k^{(I)}(R,t).
}
$$

Therefore,

$$
\boxed{
\Psi(r,R,t)
=
\sum_I\sum_k
C_k^{(I)}(t)
\chi_k^{(I)}(R,t)
\Phi_I(r;R).
}
$$

This is the central structural ansatz behind multiple-spawning methods.

---

# 22. Frozen trajectory basis functions

In FMS/AIMS, the nuclear basis functions are commonly frozen multidimensional
Gaussians, often called trajectory basis functions (TBFs).

A schematic 1D TBF is

$$
\chi_k^{(I)}(x,t)
=
N
\exp
\left[
-\frac{\alpha}{2}
(x-q_k^{(I)})^2
+
ip_k^{(I)}
(x-q_k^{(I)})
+
i\gamma_k^{(I)}
\right].
$$

The phase-space center follows approximately classical equations on electronic
surface $I$:

$$
\dot q_k^{(I)}
=
\frac{p_k^{(I)}}{m},
$$

$$
\dot p_k^{(I)}
=
-
\nabla E_I(q_k^{(I)}).
$$

But the coefficients are coupled quantum mechanically.

This distinction is crucial:

> The TBF centers can look classical while the **wavefunction represented by their
> coherent superposition is quantum mechanical.**

---

# 23. Why spawning is needed

Suppose a Gaussian on electronic state 1 approaches a nonadiabatic region.

If the exact wavefunction transfers amplitude to electronic state 2, a basis
containing no state-2 Gaussian near the same phase-space region cannot represent that
new branch efficiently.

The solution is to increase the basis dynamically:

$$
\boxed{
\chi_k^{(1)}
\longrightarrow
\chi_k^{(1)}
+
\chi_{\mathrm{new}}^{(2)}.
}
$$

This is **basis growth**, not a stochastic trajectory hop.

The parent Gaussian remains present. A child Gaussian is created on another
electronic state and the coupled coefficient equations determine how amplitude flows
between them.

---

# 24. FSSH versus spawning

The contrast is fundamental.

## FSSH

One trajectory has one active state at a time:

$$
a(t)=I.
$$

A random number decides whether the active surface changes.

Branching appears statistically across an ensemble.

## FMS/AIMS

The molecular wavefunction contains multiple coupled Gaussian basis functions.

When the basis is insufficient near a coupling region, a new basis function is
created.

Branching occurs directly in the wavefunction representation.

Therefore,

$$
\boxed{
\text{FSSH branching}
\neq
\text{FMS/AIMS spawning}.
}
$$

---

# 25. FMS and AIMS

**Full Multiple Spawning (FMS)** is the wavefunction/basis framework.

**Ab Initio Multiple Spawning (AIMS)** combines the multiple-spawning nuclear
wavefunction representation with electronic quantities evaluated on the fly from
electronic-structure theory.

Schematic loop:

```text
trajectory basis function reaches new geometry
                    |
                    v
       electronic-structure calculation
                    |
          +---------+---------+
          |         |         |
          v         v         v
        E_I(R)   gradients   couplings
          |         |         |
          +---------+---------+
                    |
                    v
       propagate TBFs and coefficients
                    |
            spawning criterion?
                    |
              yes /      \ no
                 v        v
            create TBF   continue
```

This repository implements only the **mathematical foundation and a minimal spawning
prototype**.

A complete AIMS code requires much more:

- robust electronic-structure interfaces;
- multiple electronic surfaces;
- phase/gauge-consistent couplings;
- spawning optimization;
- overlap and linear-dependence control;
- adaptive electronic/nuclear timesteps;
- energy and population bookkeeping;
- initial-condition sampling;
- convergence with respect to initial conditions and spawning thresholds.

---

# 26. Minimal two-state spawning model

For readability, the prototype uses a two-state **diabatic** Hamiltonian

$$
H_d(x)
=
\begin{pmatrix}
V_{11}(x)&V_{12}(x)\\
V_{12}(x)&V_{22}(x)
\end{pmatrix}.
$$

The state-resolved nuclear wavefunction is represented as

$$
\boldsymbol\Psi(x,t)
=
\begin{pmatrix}
\Psi_1(x,t)\\
\Psi_2(x,t)
\end{pmatrix}.
$$

Each component is expanded in Gaussian functions.

Using a diabatic model avoids derivative-coupling gauge complications while teaching
the central spawning idea:

- monitor the coupling region;
- create a child Gaussian on the other electronic state;
- solve the coupled wavefunction problem.

This should be viewed as a bridge to adiabatic FMS/AIMS, not as a substitute for it.

---

# 27. Exact split-operator reference

To judge Gaussian approximations, the repository includes a grid calculation.

For

$$
\hat H=\hat T+\hat V,
$$

the second-order Strang splitting is

$$
\boxed{
e^{-i\hat H\Delta t}
=
e^{-i\hat V\Delta t/2}
e^{-i\hat T\Delta t}
e^{-i\hat V\Delta t/2}
+
\mathcal O(\Delta t^3).
}
$$

In coordinate space,

$$
e^{-iV(x)\Delta t/2}
$$

is diagonal.

In momentum space,

$$
e^{-ip^2\Delta t/(2m)}
$$

is diagonal.

Therefore:

```text
psi(x)
  |
half potential step
  |
FFT
  |
full kinetic step
  |
inverse FFT
  |
half potential step
  |
psi(x,t+dt)
```

The global temporal error is second order:

$$
\boxed{
\text{global error}
=
\mathcal O(\Delta t^2).
}
$$

---

# 28. Hierarchy of approximations

The framework can now be read as a controlled hierarchy.

## Exact grid

No Gaussian shape restriction, but expensive in many dimensions.

## Heller TGA

One Gaussian, moving center and width.

Main approximation:

$$
V(x)
\rightarrow
V(q)+V'(q)\xi+\frac12V''(q)\xi^2.
$$

## Frozen single Gaussian

One Gaussian with fixed width.

Additional restriction: no width dynamics.

## Moving Gaussian basis

Several Gaussians; centers prescribed by trajectory equations; coefficients quantum.

## Variational multi-Gaussian

Several Gaussians; basis parameters respond to the TDSE projection.

## FMS/AIMS

Several electronic states plus adaptive Gaussian basis growth in nonadiabatic regions.

---

# 29. What convergence means for each level

## Split operator

Check

$$
\Delta t,\quad
\Delta x,\quad
x_{\max},
$$

and boundary density.

## Heller TGA

Compare to exact propagation and monitor where non-Gaussian structure develops.

## Moving Gaussian basis

Increase

$$
N_G
$$

and vary initial phase-space placement/width.

Check condition number of $S$.

## Variational multi-Gaussian

Increase $N_G$, tighten tangent derivative steps and least-squares cutoff.

Monitor TDSE residual.

## FMS/AIMS

Converge:

- number/distribution of initial conditions;
- spawning thresholds;
- basis overlap controls;
- electronic-structure accuracy;
- nuclear and electronic timesteps.

---

# 30. Scientific diagnostics used in this repository

## Norm

Grid representation:

$$
\|\Psi\|^2
=
\Delta x
\sum_n|\Psi(x_n)|^2.
$$

Moving nonorthogonal basis:

$$
\boxed{
\|\Psi\|^2=C^\dagger S C.
}
$$

## Energy

$$
E
=
\langle\Psi|\hat H|\Psi\rangle.
$$

## Wavefunction overlap with reference

Because global phase is physically irrelevant,

$$
F(t)
=
|\langle\Psi_\mathrm{ref}(t)|\Psi_\mathrm{approx}(t)\rangle|^2.
$$

## TDSE residual

For the variational dynamics,

$$
\boxed{
\epsilon_\mathrm{TDSE}
=
\left\|
i\dot\Psi-\hat H\Psi
\right\|.
}
$$

This directly measures how well the finite variational tangent space can reproduce the
exact TDSE velocity.

---

# 31. What the framework deliberately does not hide

1. A single Gaussian cannot represent arbitrary wavepacket splitting.
2. Frozen Gaussian centers can follow classical trajectories while the total
   Gaussian superposition remains quantum.
3. Nonorthogonal basis functions require overlap matrices.
4. A moving basis requires the $\tau$ term.
5. Variational parameter equations can become ill conditioned.
6. Numerical TDVP tangent derivatives are pedagogically transparent but not efficient
   for large systems.
7. A spawning demonstration is not automatically a production FMS/AIMS algorithm.
8. Adiabatic multi-state Gaussian dynamics requires careful electronic-state
   phase/gauge treatment.
9. On-the-fly AIMS requires an electronic-structure layer that is intentionally outside
   the first version of this framework.
10. Every approximation should ultimately be checked against an exact or systematically
    converged reference where one is available.

The purpose of the repository is to make each approximation visible rather than bury
it inside software machinery.
