# Detailed Derivations and Simplifications

This file collects algebra that would interrupt the flow of `THEORY.md` but is useful
when checking the implementation line by line.

---

# A. Gaussian normalization

For

$$
g(x)=N\exp[-\alpha(x-q)^2/2+ip(x-q)],
$$

the probability density is

$$
|g(x)|^2=|N|^2e^{-\alpha(x-q)^2}.
$$

Set

$$
1=|N|^2\int_{-\infty}^{\infty}e^{-\alpha y^2}dy
=|N|^2\sqrt{\frac{\pi}{\alpha}}.
$$

Therefore,

$$
\boxed{
N=\left(\frac{\alpha}{\pi}\right)^{1/4}.
}
$$

---

# B. Coordinate moments

Because the density is even around $q$,

$$
\langle x-q\rangle=0.
$$

Hence

$$
\boxed{\langle x\rangle=q.}
$$

For the variance,

$$
\langle(x-q)^2\rangle
=
\sqrt{\frac{\alpha}{\pi}}
\int y^2e^{-\alpha y^2}dy.
$$

Using

$$
\int_{-\infty}^{\infty}y^2e^{-\alpha y^2}dy
=
\frac{\sqrt\pi}{2\alpha^{3/2}},
$$

we obtain

$$
\boxed{
\langle(x-q)^2\rangle=\frac{1}{2\alpha}.
}
$$

---

# C. Momentum moments

From

$$
\partial_xg=[-\alpha(x-q)+ip]g,
$$

we get

$$
\hat p g
=
-i\partial_xg
=
[p+i\alpha(x-q)]g.
$$

Therefore,

$$
\boxed{\langle \hat p\rangle=p.}
$$

Applying $\hat p^2$ or using the second derivative gives

$$
\boxed{
\langle \hat p^2\rangle
=
p^2+\frac{\alpha}{2}.
}
$$

Thus

$$
\boxed{
\sigma_p^2=\frac{\alpha}{2}.
}
$$

Together,

$$
\sigma_x\sigma_p
=
\sqrt{\frac{1}{2\alpha}}
\sqrt{\frac{\alpha}{2}}
=
\boxed{\frac12}.
$$

---

# D. Analytic overlap of equal-width frozen Gaussians

Let

$$
g_i(x)
=
N
e^{-\alpha(x-q_i)^2/2+ip_i(x-q_i)},
$$

$$
g_j(x)
=
N
e^{-\alpha(x-q_j)^2/2+ip_j(x-q_j)}.
$$

Then

$$
S_{ij}
=
N^2
\int
e^{-\frac{\alpha}{2}[(x-q_i)^2+(x-q_j)^2]}
e^{i[p_j(x-q_j)-p_i(x-q_i)]}
dx.
$$

Define

$$
\bar q=\frac{q_i+q_j}{2},
\qquad
\Delta q=q_i-q_j,
\qquad
\Delta p=p_j-p_i.
$$

The quadratic coordinate term satisfies

$$
(x-q_i)^2+(x-q_j)^2
=
2(x-\bar q)^2+\frac{(\Delta q)^2}{2}.
$$

Thus

$$
-\frac{\alpha}{2}
[(x-q_i)^2+(x-q_j)^2]
=
-\alpha(x-\bar q)^2
-\frac{\alpha(\Delta q)^2}{4}.
$$

The momentum phase is

$$
p_j(x-q_j)-p_i(x-q_i)
=
\Delta p\,x-p_jq_j+p_iq_i.
$$

Shift

$$
y=x-\bar q.
$$

The integral becomes a standard Fourier transform of a Gaussian:

$$
\int e^{-\alpha y^2+i\Delta p\,y}dy
=
\sqrt{\frac{\pi}{\alpha}}
e^{-(\Delta p)^2/(4\alpha)}.
$$

Because

$$
N^2\sqrt{\frac{\pi}{\alpha}}=1,
$$

the result is

$$
\boxed{
S_{ij}
=
\exp
\left[
-\frac{\alpha}{4}(q_i-q_j)^2
-\frac{(p_i-p_j)^2}{4\alpha}
+\frac{i}{2}(p_i+p_j)(q_i-q_j)
\right].
}
$$

Checks:

- $S_{ii}=1$,
- $S_{ji}=S_{ij}^*$,
- $|S_{ij}|\le1$.

---

# E. Kinetic-energy action on a frozen Gaussian

Let

$$
g(x)
=
N e^{-\alpha(x-q)^2/2+ip(x-q)}.
$$

Define

$$
f(x)=-\alpha(x-q)+ip.
$$

Then

$$
g'=fg.
$$

Differentiate once more:

$$
g''
=
f'g+fg'
=
(-\alpha+f^2)g.
$$

Therefore,

$$
\boxed{
g''
=
\left(
[-\alpha(x-q)+ip]^2-\alpha
\right)g.
}
$$

The kinetic operator gives

$$
\boxed{
\hat T g
=
-\frac{1}{2m}
\left(
[-\alpha(x-q)+ip]^2-\alpha
\right)g.
}
$$

The code can therefore evaluate kinetic matrix elements either analytically through
this expression or numerically on the common grid. The framework uses the explicit
formula because it is both simple and auditable.

---

# F. Derivatives with respect to Gaussian parameters

Again use

$$
g(x;q,p,\alpha)
=
N_\alpha
e^{-\alpha(x-q)^2/2+ip(x-q)}.
$$

## F.1 Position derivative

The normalization does not depend on $q$:

$$
\boxed{
\partial_qg
=
[\alpha(x-q)-ip]g.
}
$$

## F.2 Momentum derivative

$$
\boxed{
\partial_pg
=
i(x-q)g.
}
$$

## F.3 Log-width derivative

For the variational implementation define

$$
\lambda=\ln\alpha,
\qquad
\alpha=e^\lambda.
$$

Because

$$
N_\alpha=(\alpha/\pi)^{1/4},
$$

we have

$$
\partial_\lambda\ln N_\alpha=\frac14.
$$

Also,

$$
\partial_\lambda
\left[
-\frac{\alpha}{2}(x-q)^2
\right]
=
-\frac{\alpha}{2}(x-q)^2.
$$

Therefore,

$$
\boxed{
\partial_\lambda g
=
\left[
\frac14
-
\frac{\alpha}{2}(x-q)^2
\right]g.
}
$$

This analytic result is useful for checking the numerical tangent vectors used by the
TDVP code.

---

# G. Time derivative of a frozen moving Gaussian

If $\alpha$ is fixed,

$$
\dot g
=
\partial_qg\,\dot q
+
\partial_pg\,\dot p.
$$

Substituting the derivatives above,

$$
\boxed{
\dot g
=
\left[
(\alpha(x-q)-ip)\dot q
+
i(x-q)\dot p
\right]g.
}
$$

With classical center motion,

$$
\dot q=\frac{p}{m},
\qquad
\dot p=-V'(q).
$$

This is the quantity used in

$$
\tau_{ij}
=
\langle g_i|\dot g_j\rangle.
$$

---

# H. Derivation of moving-basis norm conservation

The norm is

$$
N=C^\dagger S C.
$$

Differentiate:

$$
\dot N
=
\dot C^\dagger SC
+
C^\dagger\dot S C
+
C^\dagger S\dot C.
$$

The coefficient equation is

$$
iS\dot C=(H-i\tau)C.
$$

Hence

$$
S\dot C
=
-iHC-\tau C.
$$

Its Hermitian conjugate is

$$
\dot C^\dagger S
=
iC^\dagger H
-
C^\dagger\tau^\dagger
$$

for Hermitian $H$.

Therefore,

$$
\dot N
=
iC^\dagger HC
-
C^\dagger\tau^\dagger C
+
C^\dagger\dot S C
-iC^\dagger HC
-
C^\dagger\tau C.
$$

So

$$
\dot N
=
C^\dagger
(\dot S-\tau^\dagger-\tau)
C.
$$

But

$$
\dot S=\tau^\dagger+\tau.
$$

Thus

$$
\boxed{
\dot N=0.
}
$$

This shows precisely why the moving-basis term $\tau$ is required.

---

# I. Harmonic oscillator and exact Heller propagation

Let

$$
V(x)=\frac12m\omega^2x^2.
$$

Then

$$
V'(q)=m\omega^2q,
\qquad
V''(q)=m\omega^2.
$$

The center equations are

$$
\dot q=\frac{p}{m},
$$

$$
\dot p=-m\omega^2q.
$$

Therefore,

$$
\ddot q+\omega^2q=0,
$$

with solution

$$
\boxed{
q(t)
=
q_0\cos\omega t
+
\frac{p_0}{m\omega}\sin\omega t.
}
$$

Similarly,

$$
\boxed{
p(t)
=
p_0\cos\omega t
-
m\omega q_0\sin\omega t.
}
$$

The width equation is

$$
\dot A
=
-\frac{A^2}{m}
-
m\omega^2.
$$

For the coherent-state width

$$
A=i m\omega,
$$

we obtain

$$
-\frac{A^2}{m}-m\omega^2
=
-\frac{-m^2\omega^2}{m}
-m\omega^2
=0.
$$

Thus the coherent-state width remains constant.

This gives a useful special case:

$$
\boxed{
A_0=im\omega
\quad\Rightarrow\quad
A(t)=im\omega.
}
$$

The corresponding coordinate variance is

$$
\sigma_x^2=\frac{1}{2m\omega}.
$$

This is the harmonic-oscillator coherent-state variance.

---

# J. McLachlan TDVP derivation with real parameters

Let

$$
|\Psi\rangle=|\Psi(\theta_1,\ldots,\theta_M)\rangle
$$

with real $\theta_\mu$.

Define

$$
|D_\mu\rangle
=
\partial_{\theta_\mu}|\Psi\rangle.
$$

Then

$$
|\dot\Psi\rangle
=
\sum_\nu |D_\nu\rangle\dot\theta_\nu.
$$

The TDSE residual is

$$
|r\rangle
=
i\sum_\nu|D_\nu\rangle\dot\theta_\nu
-
H|\Psi\rangle.
$$

Minimize

$$
\mathcal R=\langle r|r\rangle
$$

with respect to each real $\dot\theta_\mu$.

Because

$$
\frac{\partial r}{\partial\dot\theta_\mu}
=
i|D_\mu\rangle,
$$

stationarity gives

$$
0
=
2\operatorname{Re}
\left\langle
iD_\mu
\middle|
r
\right\rangle.
$$

Substitute $r$:

$$
0
=
2\operatorname{Re}
\left[
\sum_\nu
\langle iD_\mu|iD_\nu\rangle
\dot\theta_\nu
-
\langle iD_\mu|H|\Psi\rangle
\right].
$$

Now

$$
\langle iD_\mu|iD_\nu\rangle
=
\langle D_\mu|D_\nu\rangle,
$$

and

$$
\operatorname{Re}
\langle iD_\mu|H|\Psi\rangle
=
\operatorname{Im}
\langle D_\mu|H|\Psi\rangle.
$$

Therefore,

$$
\boxed{
\sum_\nu
\operatorname{Re}
\langle D_\mu|D_\nu\rangle
\dot\theta_\nu
=
\operatorname{Im}
\langle D_\mu|H|\Psi\rangle.
}
$$

Or

$$
\boxed{
G\dot\theta=b.
}
$$

---

# K. Why the TDVP residual can only decrease after projection

The exact TDSE velocity is

$$
|\dot\Psi_\mathrm{exact}\rangle=-iH|\Psi\rangle.
$$

The variational manifold allows only velocities in the tangent space

$$
\mathcal T
=
\operatorname{span}
\{
|D_1\rangle,\ldots,|D_M\rangle
\}.
$$

The McLachlan solution is the least-squares projection of the exact velocity into
$\mathcal T$.

Therefore, among all possible parameter velocities,

$$
\dot\theta,
$$

the chosen one minimizes

$$
\left\|
\sum_\nu D_\nu\dot\theta_\nu
+
iH\Psi
\right\|.
$$

This gives a direct test: the optimized tangent-space residual must be no worse than
the residual obtained by setting every parameter velocity to zero.

---

# L. Two-state Gaussian wavefunction

For two electronic diabatic states,

$$
\boldsymbol\Psi(x,t)
=
\begin{pmatrix}
\Psi_1(x,t)\\
\Psi_2(x,t)
\end{pmatrix}.
$$

The Hamiltonian is

$$
\hat H
=
-\frac{1}{2m}\partial_x^2 I
+
\begin{pmatrix}
V_{11}(x)&V_{12}(x)\\
V_{12}(x)&V_{22}(x)
\end{pmatrix}.
$$

Expanding each component,

$$
\Psi_I(x,t)
=
\sum_k C_k^{(I)}g_k^{(I)}(x,t),
$$

gives

$$
\Psi
=
\sum_{I,k}
C_k^{(I)}
|g_k^{(I)}\rangle|I\rangle.
$$

The Hamiltonian matrix element between basis functions is

$$
\boxed{
H_{Ik,Jl}
=
\delta_{IJ}
\langle g_k^{(I)}|\hat T|g_l^{(J)}\rangle
+
\langle g_k^{(I)}|
V_{IJ}(x)
|g_l^{(J)}\rangle.
}
$$

For $I\ne J$, the coupling is entirely through the off-diagonal potential matrix
element in a diabatic representation.

This is the simplest algebraic setting in which to demonstrate spawning without
introducing adiabatic derivative-coupling gauge issues.
