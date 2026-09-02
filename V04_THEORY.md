# v0.4 Theory: Two-Dimensional Conical-Intersection Gaussian Dynamics

Version 0.4 introduces the first genuinely multidimensional nonadiabatic model in this
series. Its purpose is to connect four ideas that are usually learned separately:

1. multidimensional Heller Gaussian dynamics;
2. the branching-plane topology of a conical intersection;
3. vector nonadiabatic couplings and geometric phase;
4. dynamically spawned, coherently coupled Gaussian basis functions.

Atomic units are used unless otherwise stated.

---

# 1. Two-dimensional nuclear Hamiltonian

Let the nuclear coordinate be

$$
\mathbf R =
\begin{pmatrix}
x\\y
\end{pmatrix},
$$

with equal scalar nuclear mass $M$ for the two coordinates.

The diabatic Hamiltonian is

$$
\boxed{
\hat H_d
=
-\frac{1}{2M}
\left(
\frac{\partial^2}{\partial x^2}
+
\frac{\partial^2}{\partial y^2}
\right)I_2
+
V_d(x,y).
}
$$

The nuclear wavefunction is a two-component electronic spinor,

$$
\boxed{
\boldsymbol\Psi_d(\mathbf R,t)
=
\begin{pmatrix}
\psi_1(\mathbf R,t)\\
\psi_2(\mathbf R,t)
\end{pmatrix}.
}
$$

This representation is used for the exact FFT benchmark because the nuclear kinetic
operator is diagonal in electronic space.

---

# 2. Linear vibronic-coupling conical-intersection model

Use the two-state linear vibronic-coupling (LVC) potential

$$
\boxed{
V_d(x,y)
=
U(x,y)I_2
+
\begin{pmatrix}
\kappa x & \lambda y\\
\lambda y & -\kappa x
\end{pmatrix},
}
$$

where

$$
U(x,y)
=
\frac12\omega^2(x^2+y^2).
$$

Define

$$
h_z=\kappa x,
\qquad
h_x=\lambda y,
$$

and

$$
\rho
=
\sqrt{h_z^2+h_x^2}.
$$

The adiabatic energies are

$$
\boxed{
E_\pm(x,y)
=
U(x,y)\pm\rho.
}
$$

At

$$
x=y=0,
$$

we have

$$
\rho=0
$$

and therefore

$$
\boxed{
E_+(0,0)=E_-(0,0).
}
$$

The degeneracy is conical because, to first order in the branching-plane coordinates,

$$
E_+-E_-
=
2\sqrt{
\kappa^2x^2+\lambda^2y^2
}.
$$

The two coordinates have distinct physical roles:

- $x$ changes the diabatic energy difference;
- $y$ changes the diabatic interstate coupling.

This is the minimal two-coordinate topology required for an isolated two-state
conical intersection.

---

# 3. Branching-plane vectors

For a general two-state degeneracy, two first-order directions lift the degeneracy.

A diabatic model makes those directions explicit.

The energy-difference direction is

$$
\boxed{
\mathbf g
=
\nabla_{\mathbf R}
\frac{V_{11}-V_{22}}{2}.
}
$$

For the present model,

$$
\boxed{
\mathbf g=
\begin{pmatrix}
\kappa\\
0
\end{pmatrix}.
}
$$

The interstate-coupling direction is

$$
\boxed{
\mathbf h
=
\nabla_{\mathbf R}V_{12}
=
\begin{pmatrix}
0\\
\lambda
\end{pmatrix}.
}
$$

The local adiabatic gap is therefore

$$
\boxed{
\Delta E
=
2\sqrt{
(\mathbf g\cdot\mathbf R)^2
+
(\mathbf h\cdot\mathbf R)^2
}.
}
$$

This is the canonical branching-plane structure.

---

# 4. Adiabatic eigenvectors and mixing angle

Define the angle

$$
\boxed{
\theta(x,y)
=
\operatorname{atan2}(\lambda y,\kappa x).
}
$$

Then

$$
\cos\theta=\frac{\kappa x}{\rho},
\qquad
\sin\theta=\frac{\lambda y}{\rho}.
$$

One convenient real adiabatic gauge is

$$
\boxed{
|\phi_-\rangle
=
\begin{pmatrix}
-\sin(\theta/2)\\
\cos(\theta/2)
\end{pmatrix},
\qquad
|\phi_+\rangle
=
\begin{pmatrix}
\cos(\theta/2)\\
\sin(\theta/2)
\end{pmatrix}.
}
$$

Direct substitution gives

$$
V_d|\phi_\pm\rangle=E_\pm|\phi_\pm\rangle.
$$

The half-angle is crucial. If the nuclei make one full revolution around the
intersection,

$$
\theta\rightarrow\theta+2\pi.
$$

But

$$
\sin\left(\frac{\theta+2\pi}{2}\right)
=
-\sin\left(\frac{\theta}{2}\right),
$$

and similarly for the cosine. Therefore

$$
\boxed{
|\phi_\pm(\theta+2\pi)\rangle
=
-
|\phi_\pm(\theta)\rangle.
}
$$

A continuous real electronic eigenvector changes sign around a loop enclosing the CI.

---

# 5. Vector derivative coupling

Define

$$
\boxed{
\mathbf d_{-+}
=
\langle\phi_-|
\nabla_{\mathbf R}
\phi_+\rangle.
}
$$

Differentiate the upper state:

$$
\nabla|\phi_+\rangle
=
\frac12
(\nabla\theta)
\begin{pmatrix}
-\sin(\theta/2)\\
\cos(\theta/2)
\end{pmatrix}.
$$

The vector in parentheses is exactly $|\phi_-\rangle$, so

$$
\boxed{
\mathbf d_{-+}
=
\frac12\nabla\theta.
}
$$

For

$$
\theta=\operatorname{atan2}(\lambda y,\kappa x),
$$

the derivatives are

$$
\frac{\partial\theta}{\partial x}
=
-\frac{\kappa\lambda y}
{\kappa^2x^2+\lambda^2y^2},
$$

$$
\frac{\partial\theta}{\partial y}
=
\frac{\kappa\lambda x}
{\kappa^2x^2+\lambda^2y^2}.
$$

Hence

$$
\boxed{
\mathbf d_{-+}(x,y)
=
\frac{\kappa\lambda}
{2(\kappa^2x^2+\lambda^2y^2)}
\begin{pmatrix}
-y\\
x
\end{pmatrix}.
}
$$

The reverse coupling is

$$
\boxed{
\mathbf d_{+-}=-\mathbf d_{-+}.
}
$$

The $1/\rho$-type singular behavior near the CI is not a numerical accident. Individual
adiabatic electronic states cease to be uniquely defined at the degeneracy.

---

# 6. Berry/Longuet-Higgins phase from the line integral

Because

$$
\mathbf d_{-+}
=
\frac12\nabla\theta,
$$

for a closed path $\mathcal C$ that winds once around the CI,

$$
\oint_{\mathcal C}
\mathbf d_{-+}\cdot d\mathbf R
=
\frac12
\oint_{\mathcal C}
\nabla\theta\cdot d\mathbf R.
$$

The angle increases by $2\pi$, so

$$
\boxed{
\oint_{\mathcal C}
\mathbf d_{-+}\cdot d\mathbf R
=
\pi.
}
$$

This is the same topology encoded by the sign change of the real adiabatic
eigenvector.

The geometric phase is therefore not a cosmetic phase convention. It affects the
single-valuedness and interference structure of the nuclear wavefunction.

---

# 7. Parallel transport and subspace transport

For a nondegenerate real state sampled on a path, a simple discrete gauge chooses the
sign of the new eigenvector so that

$$
\boxed{
\langle\phi_i(\mathbf R_n)
|
\phi_i(\mathbf R_{n+1})
\rangle>0.
}
$$

If a closed loop encloses the CI, this continuous transport returns the final vector
with negative overlap with the initial vector.

Near an exact or numerical degeneracy, state-by-state sign alignment is not
well-defined. Instead consider an $m$-dimensional electronic subspace.

Let

$$
\Phi_{\rm ref}
$$

and

$$
\Phi_{\rm new}
$$

contain orthonormal basis vectors as columns.

Find the unitary matrix $Q$ minimizing

$$
\boxed{
\|\Phi_{\rm new}Q-\Phi_{\rm ref}\|_F.
}
$$

If

$$
\Phi_{\rm new}^\dagger\Phi_{\rm ref}
=
U\Sigma V^\dagger
$$

is an SVD, the orthogonal/unitary Procrustes solution is

$$
\boxed{
Q=UV^\dagger.
}
$$

This aligns the **subspace** without pretending that individual states are unique
inside a degenerate manifold.

---

# 8. Exact two-dimensional two-state split operator

Write

$$
\hat H_d=\hat T I_2+V_d(x,y).
$$

Second-order Strang splitting gives

$$
\boxed{
e^{-i\hat H_d\Delta t}
=
e^{-iV_d\Delta t/2}
e^{-i\hat T\Delta t}
e^{-iV_d\Delta t/2}
+
\mathcal O(\Delta t^3).
}
$$

For each grid point, diagonalize the $2\times2$ Hermitian matrix

$$
V_d(\mathbf R)=U_V\varepsilon U_V^\dagger,
$$

so

$$
e^{-iV_d\Delta t/2}
=
U_Ve^{-i\varepsilon\Delta t/2}U_V^\dagger.
$$

For the kinetic operator,

$$
\hat T
=
-\frac{1}{2M}(\partial_x^2+\partial_y^2),
$$

the momentum-space phase is

$$
\boxed{
e^{-iT\Delta t}
=
\exp\left[
-\frac{i\Delta t}{2M}
(k_x^2+k_y^2)
\right].
}
$$

A two-dimensional FFT therefore gives the exact-grid reference up to spatial and
Strang-splitting errors.

---

# 9. Multidimensional normalized Gaussian

Let the dimension be $D$. Introduce

$$
\boldsymbol\xi=\mathbf R-\mathbf q.
$$

A frozen multidimensional Gaussian with a real symmetric positive-definite width
matrix $A$ and real symmetric chirp matrix $K$ is

$$
\boxed{
g(\mathbf R)
=
N
\exp
\left[
-\frac12\boldsymbol\xi^TA\boldsymbol\xi
+
i\mathbf p^T\boldsymbol\xi
+
\frac{i}{2}\boldsymbol\xi^TK\boldsymbol\xi
\right].
}
$$

Its probability density is

$$
|g|^2
=
|N|^2
e^{-\boldsymbol\xi^TA\boldsymbol\xi}.
$$

Using the multidimensional Gaussian integral

$$
\int_{\mathbb R^D}
e^{-\mathbf x^TA\mathbf x}
d^D x
=
\frac{\pi^{D/2}}{\sqrt{\det A}},
$$

normalization gives

$$
\boxed{
N
=
\left(
\frac{\det A}{\pi^D}
\right)^{1/4}.
}
$$

The coordinate covariance matrix is

$$
\boxed{
\operatorname{Cov}(\mathbf R)
=
\frac12A^{-1}.
}
$$

---

# 10. Multidimensional Gaussian kinetic energy

Define the complex symmetric matrix

$$
\boxed{
Z=A-iK.
}
$$

Then

$$
g
=
N
\exp
\left[
-\frac12\xi^TZ\xi+i\mathbf p^T\xi
\right].
$$

The gradient is

$$
\boxed{
\nabla g
=
(-Z\xi+i\mathbf p)g.
}
$$

Define

$$
\mathbf f=-Z\xi+i\mathbf p.
$$

For each coordinate,

$$
\partial_\alpha^2g
=
(f_\alpha^2-Z_{\alpha\alpha})g.
$$

Summing,

$$
\boxed{
\nabla^2g
=
\left(
\mathbf f^T\mathbf f-\operatorname{Tr}Z
\right)g.
}
$$

Therefore,

$$
\boxed{
\hat T g
=
-\frac{1}{2M}
\left(
\mathbf f^T\mathbf f-\operatorname{Tr}Z
\right)g.
}
$$

This exact expression is used by the 2D Gaussian matrix-element code.

---

# 11. Equal-width multidimensional Gaussian overlap

For two zero-chirp Gaussians sharing the same real positive-definite width matrix $A$,

$$
g_i=g(\mathbf q_i,\mathbf p_i,A),
\qquad
g_j=g(\mathbf q_j,\mathbf p_j,A),
$$

define

$$
\Delta\mathbf q=\mathbf q_i-\mathbf q_j,
\qquad
\Delta\mathbf p=\mathbf p_i-\mathbf p_j.
$$

Completing the square gives

$$
\boxed{
S_{ij}
=
\exp
\left[
-\frac14\Delta q^TA\Delta q
-\frac14\Delta p^TA^{-1}\Delta p
+\frac{i}{2}
(\mathbf p_i+\mathbf p_j)^T
(\mathbf q_i-\mathbf q_j)
\right].
}
$$

The magnitude is controlled by phase-space separation,

$$
\boxed{
|S_{ij}|
=
\exp
\left[
-\frac14\Delta q^TA\Delta q
-\frac14\Delta p^TA^{-1}\Delta p
\right].
}
$$

This naturally supplies an overlap-based spawning redundancy criterion.

---

# 12. Multidimensional Heller thawed Gaussian

Use the Heller phase convention

$$
\boxed{
\Psi(\mathbf R,t)
=
\exp
\left\{
i\left[
\frac12\xi^TA_t\xi
+
\mathbf p_t^T\xi
+
\gamma_t
\right]
\right\},
}
$$

where now $A_t$ is a **complex symmetric matrix** with positive-definite imaginary
part.

Expand the potential locally:

$$
V(\mathbf R)
\approx
V(\mathbf q)
+
\nabla V(\mathbf q)^T\xi
+
\frac12\xi^T
H_V(\mathbf q)
\xi,
$$

where

$$
H_V=\nabla\nabla V
$$

is the Hessian.

Matching quadratic, linear, and constant powers of $\xi$ gives

$$
\boxed{
\dot{\mathbf q}
=
\frac{\mathbf p}{M},
}
$$

$$
\boxed{
\dot{\mathbf p}
=
-\nabla V(\mathbf q),
}
$$

$$
\boxed{
\dot A
=
-\frac{1}{M}A^2-H_V(\mathbf q),
}
$$

and

$$
\boxed{
\dot\gamma
=
\frac{\mathbf p^T\mathbf p}{2M}
-
V(\mathbf q)
+
\frac{i}{2M}\operatorname{Tr}A.
}
$$

The matrix Riccati equation allows widths along different coordinates to breathe and
become correlated.

---

# 13. Symmetry of the width matrix

If

$$
A(0)=A(0)^T
$$

and the potential Hessian is symmetric,

$$
H_V=H_V^T,
$$

then

$$
\dot A^T
=
-\frac{1}{M}(A^T)^2-H_V^T
=
\dot A.
$$

Thus a symmetric initial width remains symmetric under the exact TGA equations.

The code symmetrizes only at roundoff level after each numerical step; that operation
is a numerical safeguard, not a change to the underlying equations.

---

# 14. Moving Gaussian basis for multiple electronic states

Let a trajectory basis function (TBF) carry:

- electronic state $a_k$;
- center $\mathbf q_k$;
- momentum $\mathbf p_k$;
- width matrix $A_k$.

The molecular wavefunction is approximated by

$$
\boxed{
|\Psi(t)\rangle
=
\sum_{k=1}^{N_{\rm TBF}}
C_k(t)|G_k(t)\rangle.
}
$$

The TBFs are nonorthogonal. Define

$$
S_{ij}
=
\langle G_i|G_j\rangle,
$$

$$
H_{ij}
=
\langle G_i|\hat H|G_j\rangle,
$$

$$
T_{ij}
=
\langle G_i|\dot G_j\rangle.
$$

Projection of the TDSE gives

$$
\boxed{
iS\dot C
=
(H-iT)C.
}
$$

This is the same moving-basis equation derived in v0.1, now applied to a dynamically
growing multidimensional nonadiabatic basis.

---

# 15. Why coefficient propagation matters

Spawning alone is not quantum dynamics.

If a child Gaussian is created but its coefficient is never coupled to the parent,
the basis merely grows geometrically.

In v0.4, after a child is created with initial coefficient

$$
C_{\rm child}=0,
$$

the enlarged matrices $S,H,T$ are rebuilt and the coupled coefficient equation

$$
iS\dot C=(H-iT)C
$$

is propagated.

Therefore population transfer occurs through Hamiltonian matrix elements, not by
manually assigning a population to the child.

This is the essential difference between a spawning **demonstration** and a coupled
multiple-spawning propagation.

---

# 16. Adiabatic trajectory guidance

For a TBF associated with adiabatic state $a$,

$$
\boxed{
\dot{\mathbf q}_k
=
\frac{\mathbf p_k}{M},
}
$$

$$
\boxed{
\dot{\mathbf p}_k
=
-\nabla E_a(\mathbf q_k).
}
$$

For the LVC model,

$$
\nabla E_\pm
=
\omega^2
\begin{pmatrix}
x\\y
\end{pmatrix}
\pm
\frac{1}{\rho}
\begin{pmatrix}
\kappa^2x\\
\lambda^2y
\end{pmatrix},
$$

away from the CI.

This is the trajectory-guidance approximation used in the adiabatic spawned-basis
prototype.

---

# 17. Vector-NAC spawning indicator

For parent state $a$ and target $b$, define

$$
\boxed{
\eta_{ab}
=
\left|
\dot{\mathbf q}\cdot
\mathbf d_{ab}(\mathbf q)
\right|.
}
$$

Since

$$
\dot{\mathbf q}=\frac{\mathbf p}{M},
$$

$$
\boxed{
\eta_{ab}
=
\left|
\frac{\mathbf p}{M}\cdot
\mathbf d_{ab}(\mathbf q)
\right|.
}
$$

A spawn is considered when

$$
\eta_{ab}>\eta_{\rm spawn}.
$$

The implementation additionally blocks a spawn if an existing target-state Gaussian
already has large phase-space overlap with the proposed child.

---

# 18. Energy-conserving child momentum along the NAC

Let

$$
\hat{\mathbf n}
=
\frac{\mathbf d_{ab}}
{|\mathbf d_{ab}|}.
$$

Place the child at the parent position and write

$$
\mathbf p_b
=
\mathbf p_a+\eta\hat{\mathbf n}.
$$

Require local energy conservation:

$$
\frac{|\mathbf p_a+\eta\hat n|^2}{2M}
+
E_b
=
\frac{|\mathbf p_a|^2}{2M}
+
E_a.
$$

This gives

$$
\eta^2
+
2(\mathbf p_a\cdot\hat n)\eta
+
2M(E_b-E_a)
=
0.
$$

Therefore,

$$
\boxed{
\eta
=
-(\mathbf p_a\cdot\hat n)
\pm
\sqrt{
(\mathbf p_a\cdot\hat n)^2
-
2M(E_b-E_a)
}.
}
$$

v0.4 chooses the real root with the smallest momentum change $|\eta|$.

If the discriminant is negative, this specific local energy-conserving initialization
has no real solution and the spawn is rejected.

This is a transparent prototype rule, not the full optimal-spawning procedure used in
production AIMS.

---

# 19. Basis growth without discontinuously changing the wavefunction

Suppose the old basis is

$$
\{G_1,\ldots,G_N\}
$$

with coefficient vector

$$
C=
(C_1,\ldots,C_N)^T.
$$

After spawning $G_{N+1}$, initialize

$$
\boxed{
C'
=
(C_1,\ldots,C_N,0)^T.
}
$$

Then immediately after insertion,

$$
\sum_{i=1}^{N+1}C_i'G_i
=
\sum_{i=1}^{N}C_iG_i.
$$

So basis growth alone does not alter the represented wavefunction.

Only subsequent coupled propagation transfers amplitude into the child.

---

# 20. Numerical stability of the nonorthogonal basis

The overlap matrix may become ill-conditioned when two Gaussians are nearly redundant.

Monitor

$$
\boxed{
\kappa(S)
=
\frac{\sigma_{\max}(S)}
{\sigma_{\min}(S)}.
}
$$

v0.4 uses:

- overlap blocking before spawning;
- linear solves rather than explicit $S^{-1}$;
- a condition-number diagnostic.

A production multiple-spawning code requires more sophisticated basis pruning and
regularization.

---

# 21. What v0.4 is and is not

v0.4 **does implement**:

- a two-dimensional CI model;
- exact two-state 2D FFT dynamics;
- matrix-width multidimensional Gaussians;
- multidimensional Heller TGA;
- analytic vector derivative couplings;
- Berry-phase line-integral and sign-change diagnostics;
- unitary/orthogonal subspace alignment;
- adiabatic TBF trajectory guidance;
- vector-NAC spawning;
- dynamic basis growth;
- genuinely coupled nonorthogonal coefficient propagation.

v0.4 **does not claim** to be a production AIMS program.

Production AIMS/FMS additionally requires, among other things:

- ab initio multidimensional electronic structure;
- robust spawning optimization;
- mature saddle-point/independent-first-generation approximations;
- adaptive electronic and nuclear timesteps;
- initial-condition ensembles;
- large-scale basis management;
- full convergence studies.

The purpose of v0.4 is to make the mathematical bridge to those methods explicit and
testable before adding that engineering complexity.
