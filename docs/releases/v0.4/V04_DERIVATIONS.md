# v0.4 Detailed Derivations

This file expands algebra used in `V04_THEORY.md`.

---

## A. Adiabatic energies of the LVC Hamiltonian

Ignoring the common scalar $U$, define

$$
W=
\begin{pmatrix}
\kappa x & \lambda y\\
\lambda y & -\kappa x
\end{pmatrix}.
$$

The characteristic polynomial is

$$
\det(W-\epsilon I)
=
(\kappa x-\epsilon)(-\kappa x-\epsilon)-\lambda^2y^2.
$$

Thus

$$
\epsilon^2-\kappa^2x^2-\lambda^2y^2=0,
$$

so

$$
\epsilon_\pm
=
\pm\sqrt{\kappa^2x^2+\lambda^2y^2}.
$$

Adding $U$ gives

$$
\boxed{
E_\pm=U\pm\rho.
}
$$

---

## B. Gradient of the adiabatic energies

For

$$
\rho=(\kappa^2x^2+\lambda^2y^2)^{1/2},
$$

$$
\partial_x\rho
=
\frac{\kappa^2x}{\rho},
$$

$$
\partial_y\rho
=
\frac{\lambda^2y}{\rho}.
$$

Since

$$
\nabla U
=
\omega^2(x,y)^T,
$$

$$
\boxed{
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
\end{pmatrix}.
}
$$

---

## C. Derivative coupling from the mixing angle

Let

$$
|\phi_+\rangle=
\begin{pmatrix}
c\\s
\end{pmatrix},
\qquad
|\phi_-\rangle=
\begin{pmatrix}
-s\\c
\end{pmatrix},
$$

where

$$
c=\cos(\theta/2),
\qquad
s=\sin(\theta/2).
$$

Then

$$
\nabla|\phi_+\rangle
=
\frac12(\nabla\theta)
\begin{pmatrix}
-s\\c
\end{pmatrix}
=
\frac12(\nabla\theta)|\phi_-\rangle.
$$

Therefore

$$
\boxed{
\langle\phi_-|\nabla\phi_+\rangle
=
\frac12\nabla\theta.
}
$$

For $\theta=\operatorname{atan2}(Y,X)$,

$$
d\theta
=
\frac{X\,dY-Y\,dX}{X^2+Y^2}.
$$

Take

$$
X=\kappa x,
\qquad
Y=\lambda y.
$$

Then

$$
dX=\kappa dx,
\qquad
dY=\lambda dy,
$$

so

$$
d\theta
=
\frac{
\kappa\lambda x\,dy
-
\kappa\lambda y\,dx
}{
\kappa^2x^2+\lambda^2y^2
}.
$$

Hence

$$
\boxed{
\nabla\theta
=
\frac{\kappa\lambda}
{\kappa^2x^2+\lambda^2y^2}
\begin{pmatrix}
-y\\x
\end{pmatrix}.
}
$$

---

## D. Geometric-phase line integral on a circle

For the isotropic case

$$
\kappa=\lambda,
$$

write

$$
x=r\cos\varphi,
\qquad
y=r\sin\varphi.
$$

Then

$$
\mathbf d_{-+}
=
\frac{1}{2r}
\hat{\boldsymbol\varphi}.
$$

Also

$$
d\mathbf R=r\,d\varphi\,\hat{\boldsymbol\varphi}.
$$

Therefore

$$
\oint\mathbf d_{-+}\cdot d\mathbf R
=
\int_0^{2\pi}
\frac{1}{2r}r\,d\varphi
=
\boxed{\pi}.
$$

The result depends on winding number, not radius.

---

## E. Multidimensional Gaussian normalization

For

$$
g(\mathbf R)
=
N
e^{-\frac12\xi^TA\xi+i(\cdots)},
$$

$$
|g|^2
=
|N|^2e^{-\xi^TA\xi}.
$$

The Gaussian integral is

$$
\int e^{-\xi^TA\xi}d^D\xi
=
\frac{\pi^{D/2}}{\sqrt{\det A}}.
$$

Set the norm equal to one:

$$
1
=
|N|^2
\frac{\pi^{D/2}}{\sqrt{\det A}}.
$$

Therefore,

$$
|N|^2
=
\frac{\sqrt{\det A}}{\pi^{D/2}},
$$

and

$$
\boxed{
N
=
\left(
\frac{\det A}{\pi^D}
\right)^{1/4}.
}
$$

---

## F. Multidimensional kinetic operator

Let

$$
Z=A-iK,
$$

and

$$
g
=
N
e^{-\frac12\xi^TZ\xi+i p^T\xi}.
$$

The derivative with respect to coordinate $\alpha$ is

$$
\partial_\alpha g
=
f_\alpha g,
$$

where

$$
f_\alpha
=
-(Z\xi)_\alpha+ip_\alpha.
$$

Differentiate once more:

$$
\partial_\alpha^2g
=
(\partial_\alpha f_\alpha)g
+
f_\alpha\partial_\alpha g.
$$

Since

$$
\partial_\alpha f_\alpha=-Z_{\alpha\alpha},
$$

$$
\partial_\alpha^2g
=
(f_\alpha^2-Z_{\alpha\alpha})g.
$$

Summing over coordinates,

$$
\boxed{
\nabla^2g
=
(f^Tf-\operatorname{Tr}Z)g.
}
$$

---

## G. Multidimensional Heller width equation

Use

$$
S(\mathbf R,t)
=
\frac12\xi^TA\xi+p^T\xi+\gamma.
$$

For the kinetic term,

$$
\frac{1}{2M}(\nabla S)^T(\nabla S)
=
\frac{1}{2M}
(A\xi+p)^T(A\xi+p).
$$

The quadratic coefficient is

$$
\frac{1}{2M}\xi^TA^2\xi.
$$

The local potential contributes

$$
\frac12\xi^TH_V\xi.
$$

The time derivative contributes

$$
-\frac12\xi^T\dot A\xi.
$$

Equating quadratic forms gives

$$
-\dot A
=
\frac{1}{M}A^2+H_V,
$$

so

$$
\boxed{
\dot A=-M^{-1}A^2-H_V.
}
$$

For a diagonal scalar mass $M$, $M^{-1}$ is simply $1/M$.

---

## H. NAC-direction energy-conserving spawn

Take

$$
p_b=p_a+\eta n,
\qquad
|n|=1.
$$

Energy conservation requires

$$
\frac{(p_a+\eta n)^T(p_a+\eta n)}{2M}+E_b
=
\frac{p_a^Tp_a}{2M}+E_a.
$$

Expand:

$$
p_a^Tp_a+2\eta p_a\cdot n+\eta^2
+
2ME_b
=
p_a^Tp_a+2ME_a.
$$

Cancel the common kinetic term:

$$
\eta^2
+
2(p_a\cdot n)\eta
+
2M(E_b-E_a)
=
0.
$$

The discriminant is

$$
\boxed{
D
=
(p_a\cdot n)^2-2M(E_b-E_a).
}
$$

Real roots exist only for $D\ge0$.

---

## I. Continuity under basis insertion

Old wavefunction:

$$
\Psi_{\mathrm{old}}
=
\sum_{i=1}^NC_iG_i.
$$

New basis:

$$
\{G_1,\ldots,G_N,G_{N+1}\}.
$$

Choose

$$
C_{N+1}=0.
$$

Then

$$
\Psi_{\mathrm{new}}
=
\sum_{i=1}^NC_iG_i+0\cdot G_{N+1}
=
\boxed{\Psi_{\mathrm{old}}}.
$$

Thus spawning changes the variational space but not the instantaneous represented
state.
