# v0.25.1 multi-Gaussian TDVP: mathematical specification

## 1. Released Hamiltonian

The nuclear coordinate is `x`, its constant mass is `m>0`, and the electronic
space has `N_s` complete spin states in one coordinate-independent orthonormal
frame. The Hamiltonian is

$$
\hat H=-\frac{1}{2m}\frac{d^2}{dx^2}\mathbf 1_{N_s}
       +H_0+xH_1+x^2H_2,
\qquad H_k=H_k^\dagger.
$$

The quadratic restriction makes the released matrix-element oracle analytic and
auditable. Provider intake verifies the Hamiltonian and physical derivative on an
independent coordinate stencil, checks zero electronic connection, and requires a
complete multiplet provenance record.

## 2. Frozen-width spinor ansatz

For packet `I`,

$$
g_I(x)=\left(\frac{\alpha_I}{\pi}\right)^{1/4}
\exp\left[-\frac{\alpha_I}{2}(x-q_I)^2
          +ip_I(x-q_I)\right],\qquad \alpha_I>0.
$$

The total wavefunction is

$$
|\Psi(\theta)\rangle=\sum_{I=1}^{N_g}g_I(x)
\sum_{a=1}^{N_s}C_{Ia}|a\rangle.
$$

The real parameter vector is stored in the exact order

$$
\theta=(\operatorname{Re}C_{11},\ldots,\operatorname{Re}C_{N_gN_s},
        \operatorname{Im}C_{11},\ldots,\operatorname{Im}C_{N_gN_s},
        q_1,\ldots,q_{N_g},p_1,\ldots,p_{N_g}).
$$

Thus the metric dimension is

$$P=2N_gN_s+2N_g.$$

Widths are immutable model parameters in this release: `alpha_dot=0`.

## 3. Tangent vectors

The coefficient tangents are

$$
\frac{\partial\Psi}{\partial\operatorname{Re}C_{Ia}}=g_I|a\rangle,
\qquad
\frac{\partial\Psi}{\partial\operatorname{Im}C_{Ia}}=i g_I|a\rangle.
$$

Writing `y_I=x-q_I`, the packet tangents are

$$
\frac{\partial g_I}{\partial q_I}
=[\alpha_I y_I-ip_I]g_I,
\qquad
\frac{\partial g_I}{\partial p_I}=iy_Ig_I.
$$

Consequently every tangent is a degree-zero or degree-one polynomial times one
Gaussian and one electronic vector.

## 4. McLachlan equations

McLachlan variation minimizes the instantaneous Schrödinger defect

$$
\mathcal F(\dot\theta)=
\|\dot\Psi+i\hat H\Psi\|^2
$$

over real `theta_dot`. Differentiating with respect to each velocity gives

$$
\sum_\nu G_{\mu\nu}\dot\theta_\nu=b_\mu,
$$

where

$$
G_{\mu\nu}=\operatorname{Re}
\langle\partial_\mu\Psi|\partial_\nu\Psi\rangle,
\qquad
b_\mu=\operatorname{Im}
\langle\partial_\mu\Psi|\hat H\Psi\rangle.
$$

`G` is real, symmetric, and positive semidefinite. It can be singular because
different parameter variations can represent the same physical tangent—for example,
two exactly coincident packets.

## 5. Exact Gaussian moments

For a bra/ket packet pair, define

$$
M_n^{IJ}=\langle g_I|x^n|g_J\rangle,
\qquad S_{IJ}=M_0^{IJ}.
$$

Completing the square gives a generally complex centroid `mu_IJ` and covariance
`v_IJ`. The moments required by v0.25.1 are

$$
M_0=S,
\quad M_1=S\mu,
\quad M_2=S(\mu^2+v),
\quad M_3=S(\mu^3+3\mu v).
$$

The kinetic operator on ket packet `J` is

$$
-\frac{1}{2m}g_J''=
\left[\frac{p_J^2+\alpha_J}{2m}
+\frac{i\alpha_Jp_J}{m}y_J
-\frac{\alpha_J^2}{2m}y_J^2\right]g_J.
$$

Therefore all overlap/Hamiltonian terms require moments through degree two, while a
degree-one bra tangent times the degree-two Hamiltonian action requires degree three.
No spatial grid or finite difference enters the production metric.

## 6. Norm and energy

The nonorthogonal generalized norm and variational energy are

$$
N=C^\dagger S C,
\qquad
E=\frac{C^\dagger H C}{C^\dagger S C},
$$

where combined indices `(I,a)` flatten packet and electronic labels. Initial states
are explicitly normalized. Every step records and independently recomputes both
endpoint quantities.

## 7. SVD metric solve

For

$$G=U\,\operatorname{diag}(\sigma_i)V^T,$$

v0.25.1 uses

$$
\tau=\max(\tau_{\rm abs},\tau_{\rm rel}\sigma_{\max}),
\qquad
\sigma_i^+=\begin{cases}1/\sigma_i,&\sigma_i>\tau,\\0,&\sigma_i\le\tau.
\end{cases}
$$

Before accepting a rank-deficient solve, the discarded left-singular subspace must
be dynamically compatible:

$$
\frac{\|U_0^Tb\|_2}{\max(\|b\|_2,10^{-30})}
\le\epsilon_{\rm null}.
$$

The minimum-norm velocity is

$$\dot\theta=V\,\operatorname{diag}(\sigma_i^+)U^Tb.$$

The retained condition number and relative residual

$$
\kappa=\sigma_{\max}/\sigma_{\min,\rm retained},
\qquad
r_G=\|G\dot\theta-b\|_2/\max(\|b\|_2,10^{-30})
$$

must also pass their release gates.

## 8. Fully implicit midpoint

Let `v(theta)` denote the SVD-resolved TDVP velocity. A signed step `h` solves

$$
R(\theta_{n+1})=\theta_{n+1}-\theta_n
-h\,v\left(\frac{\theta_n+\theta_{n+1}}{2}\right)=0.
$$

The initial nonlinear guess is the explicit tangent predictor

$$\theta_{n+1}^{(0)}=\theta_n+h\,v(\theta_n).$$

The production solve uses `scipy.optimize.root(method="hybr")`. Acceptance requires
both solver success and an independently recomputed residual norm below tolerance.
Using `-h` with the endpoint as the new start yields the validated signed reversal.

## 9. Covariances

Relabeling packets by a permutation must permute `q`, `p`, widths, and coefficient
rows without changing the physical trajectory. For a constant electronic unitary
`U`,

$$
H_k'=U^\dagger H_kU,
\qquad C'=CU^*,
$$

which leaves the represented spinor wavefunction invariant. Both covariances are
validated to nonlinear-solver precision. Coordinate-dependent `U(x)` would introduce
connection and kinetic terms and is intentionally outside this fixed-frame release.

## 10. Exact boundary

The phrases “multi-Gaussian TDVP” and “full variational metric layer” in v0.25.1
mean the released frozen-width, one-dimensional, fixed-frame, quadratic contract.
They do not mean adaptive widths, multidimensional vMCG, AIMS spawning/pruning, or
trajectory-ready molecular electronic structure.

