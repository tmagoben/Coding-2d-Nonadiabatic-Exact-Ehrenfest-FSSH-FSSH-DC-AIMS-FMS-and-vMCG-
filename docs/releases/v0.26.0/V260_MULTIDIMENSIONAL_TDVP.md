# v0.26.0 Multidimensional Gaussian TDVP

## 1. Frozen conventions

Atomic units are used.  Nuclear coordinates are column vectors
`R in R^D`.  The fixed global diabatic electronic frame contains `S` states.  The
Hamiltonian is

$$
\hat H=-\frac12\nabla^T M^{-1}\nabla\,\mathbf 1_S+\mathbf V(\mathbf R),
$$

where `M` is real symmetric positive definite and

$$
\mathbf V(\mathbf R)=\mathbf H_0+
\sum_a \mathbf H_{1,a}R_a+
\sum_{ab}\mathbf H_{2,ab}R_aR_b.
$$

There is no hidden factor of one half in `H2`; `H2[a,b]=H2[b,a]`.  Every
electronic coefficient matrix is Hermitian.

## 2. Ansatz

The complete spinor is

$$
\Psi_s(\mathbf R,t)=\sum_{I=1}^{G}C_{Is}(t)g_I(\mathbf R,t),
$$

with

$$
g_I(\mathbf R)=\prod_{a=1}^{D}
\left(\frac{\alpha_{Ia}}{\pi}\right)^{1/4}
\exp\left[-\frac{z_{Ia}}2(R_a-q_{Ia})^2
+ip_{Ia}(R_a-q_{Ia})\right],
$$

$$
z_{Ia}=\alpha_{Ia}-i\beta_{Ia},
\qquad \alpha_{Ia}=e^{\eta_{Ia}}>0.
$$

The real parameter order is

$$
\theta=(\operatorname{Re}C,\operatorname{Im}C,
q,p,\eta,\beta),
$$

with packet-major ordering inside each block.  Its dimension is

$$
P=2GS+4GD.
$$

## 3. Exact cross moments

For a bra packet `I` and ket packet `J`, define coordinatewise

$$
Z_a=z_{Ia}^*+z_{Ja},
\qquad
L_a=z_{Ia}^*q_{Ia}+z_{Ja}q_{Ja}+i(p_{Ja}-p_{Ia}).
$$

The product is a complex normal distribution with

$$
\mu_a=\frac{L_a}{Z_a},
\qquad
\Sigma_a=\frac1{Z_a}.
$$

The overlap is

$$
S_{IJ}=\exp\sum_a\left[
\frac14\log\alpha_{Ia}+\frac14\log\alpha_{Ja}
+\frac12\log2-\frac12\log Z_a
-\frac12z_{Ia}^*q_{Ia}^2-\frac12z_{Ja}q_{Ja}^2
+ip_{Ia}q_{Ia}-ip_{Ja}q_{Ja}
+\frac{L_a^2}{2Z_a}
\right].
$$

Raw one-coordinate moments are

$$
\begin{aligned}
m_0&=1,\\
m_1&=\mu,\\
m_2&=\mu^2+\Sigma,\\
m_3&=\mu^3+3\mu\Sigma,\\
m_4&=\mu^4+6\mu^2\Sigma+3\Sigma^2.
\end{aligned}
$$

Because the released width matrices are diagonal, a multivariate monomial moment
factorizes into `S_IJ` times the product of these raw moments.

## 4. Kinetic polynomial

Let

$$
f_a(\mathbf R)=-z_{Ja}(R_a-q_{Ja})+ip_{Ja}.
$$

Then

$$
\partial_a\partial_b g_J=
\left[f_af_b-\delta_{ab}z_{Ja}\right]g_J,
$$

and therefore

$$
\hat T g_J=-\frac12\sum_{ab}(M^{-1})_{ab}
\left[f_af_b-\delta_{ab}z_{Ja}\right]g_J.
$$

This is a degree-two polynomial times `g_J`.  Multiplication by the quadratic
potential remains degree two, while contraction with a degree-two shape tangent
requires moments only through degree four.

## 5. Tangents

For each packet and coordinate,

$$
\partial_{q_a}g=[z_a(R_a-q_a)-ip_a]g,
$$

$$
\partial_{p_a}g=i(R_a-q_a)g,
$$

$$
\partial_{\eta_a}g=
\left[\frac14-\frac{\alpha_a}{2}(R_a-q_a)^2\right]g,
$$

$$
\partial_{\beta_a}g=\frac{i}{2}(R_a-q_a)^2g.
$$

Coefficient real and imaginary tangents are `g|s>` and `i g|s>`.  Shape tangents
are multiplied by the complete electronic coefficient row `C[I,:]`.

## 6. McLachlan system

For real parameters,

$$
G_{\mu\nu}=\operatorname{Re}
\langle\partial_\mu\Psi|\partial_\nu\Psi\rangle,
\qquad
b_\mu=\operatorname{Im}
\langle\partial_\mu\Psi|\hat H\Psi\rangle,
$$

$$
G\dot\theta=b.
$$

A full SVD is used:

$$
G=U\Sigma V^T,
\qquad
\dot\theta=V\Sigma^+U^Tb.
$$

The retained rank, nullity, cutoff, condition number, projected-null right-hand
side, linear residual, and velocity norm are recorded.  An incompatible null-space
component fails closed; diagonal loading is not used.

## 7. Implicit midpoint

One signed step solves

$$
F(\theta_{n+1})=\theta_{n+1}-\theta_n
-\Delta t\,f\!\left(\frac{\theta_n+\theta_{n+1}}2\right)=0.
$$

The explicit TDVP predictor initializes the nonlinear solve.  The stored receipt
contains the final nonlinear residual, function evaluations, midpoint metric/SVD,
endpoint norm and energy, and maximum logarithmic-width change.

## 8. Exact one-dimensional reduction

For `D=1`, the new cross moments, overlap, Hamiltonian, metric, right-hand side, and
velocity reproduce v0.25.2 numerically.  The validation maxima are respectively
approximately

- overlap: `6.94e-18`;
- Hamiltonian: `1.74e-18`;
- metric: `9.70e-17`;
- right-hand side: `2.60e-18`;
- velocity: `1.47e-16`.

This is a direct numerical reduction, not a qualitative similarity test.

## 9. Covariance boundary

Packet permutation, constant electronic `U(S)` gauge transformations, and signed
coordinate permutations are validated.  A general rotation of an anisotropic
diagonal width produces off-diagonal width elements, so arbitrary rotational
covariance is outside this manifold.  Full complex symmetric widths are required
for that next step.
