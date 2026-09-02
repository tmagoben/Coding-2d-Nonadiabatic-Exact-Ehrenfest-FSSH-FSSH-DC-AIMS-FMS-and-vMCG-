# v0.10 Detailed Derivations

---

## A. Wigner function of the frozen multidimensional Gaussian

Take

$$
\psi(q)
=
N
\exp\left[
-\frac12(q-q_0)^TA(q-q_0)
+
ip_0^T(q-q_0)
\right],
$$

with real symmetric positive-definite $A$.

The Wigner transform is

$$
W(q,p)
=
\frac{1}{(2\pi)^D}
\int
d\xi\,
e^{-ip^T\xi}
\psi^*(q-\xi/2)
\psi(q+\xi/2).
$$

The product of wavefunctions gives

$$
\psi^*(q-\xi/2)\psi(q+\xi/2)
\propto
\exp[
-(q-q_0)^TA(q-q_0)
-\frac14\xi^TA\xi
+
ip_0^T\xi
].
$$

Therefore the Fourier integral is Gaussian:

$$
\int d\xi\,
\exp[
-\frac14\xi^TA\xi
-i(p-p_0)^T\xi
]
\propto
\exp[
-(p-p_0)^TA^{-1}(p-p_0)
].
$$

Hence

$$
\boxed{
W(q,p)
\propto
\exp[
-(q-q_0)^TA(q-q_0)
-(p-p_0)^TA^{-1}(p-p_0)
].
}
$$

Comparing with the multivariate normal form

$$
e^{-\frac12x^T\Sigma^{-1}x}
$$

gives

$$
\boxed{
\Sigma_q=\frac12A^{-1},
\qquad
\Sigma_p=\frac12A.
}
$$

---

## B. Graph reduced electronic density matrix

Write

$$
|\Psi(q)\rangle
=
\sum_iC_i g_i(q)|v_i\rangle,
$$

where all $|v_i\rangle$ have been transported to the same electronic frame.

Then

$$
|\Psi\rangle\langle\Psi|
=
\sum_{ij}
C_iC_j^*
g_i(q)g_j^*(q)
|v_i\rangle\langle v_j|.
$$

Integrate over $q$:

$$
\rho_e
=
\sum_{ij}
C_iC_j^*
\left[
\int g_i(q)g_j^*(q)dq
\right]
|v_i\rangle\langle v_j|.
$$

But

$$
\int g_i g_j^*dq
=
\langle g_j|g_i\rangle.
$$

Therefore

$$
\boxed{
\rho_e
=
\sum_{ij}
C_iC_j^*
S^{\rm nuc}_{ji}
|v_i\rangle\langle v_j|.
}
$$

---

## C. Trace equals generalized norm

Take the trace:

$$
\operatorname{Tr}\rho_e
=
\sum_{ij}
C_iC_j^*
S^{\rm nuc}_{ji}
\langle v_j|v_i\rangle.
$$

Relabel $i\leftrightarrow j$:

$$
\operatorname{Tr}\rho_e
=
\sum_{ij}
C_i^*C_j
S^{\rm nuc}_{ij}
\langle v_i|v_j\rangle.
$$

The full graph-Gaussian overlap is

$$
S_{ij}
=
S^{\rm nuc}_{ij}
\langle v_i|v_j\rangle.
$$

Hence

$$
\boxed{
\operatorname{Tr}\rho_e
=
C^\dagger SC.
}
$$

This identity is the fundamental normalization test for the new observable.

---

## D. Exact fixed-frame reduced density

For the exact diabatic wavefunction,

$$
|\Psi_d(R)\rangle,
$$

define

$$
\rho_d
=
\int
|\Psi_d(R)\rangle
\langle\Psi_d(R)|
dR.
$$

Let a fixed unitary electronic frame be $U_r$ with

$$
|\Psi_d\rangle
=
U_r|\Psi_r\rangle.
$$

Then

$$
|\Psi_r\rangle
=
U_r^\dagger|\Psi_d\rangle.
$$

Therefore

$$
\rho_r
=
\int
U_r^\dagger
|\Psi_d\rangle\langle\Psi_d|
U_r\,dR.
$$

Because $U_r$ is constant with respect to $R$ in this fixed-frame observable,

$$
\boxed{
\rho_r
=
U_r^\dagger\rho_dU_r.
}
$$

---

## E. Purity and Schmidt structure

For a normalized bipartite pure state,

$$
|\Psi\rangle
=
\sum_k
\sqrt{\lambda_k}
|u_k\rangle_{\rm nuc}
|v_k\rangle_{\rm el},
$$

the electronic reduced density matrix is

$$
\rho_e
=
\sum_k
\lambda_k
|v_k\rangle\langle v_k|.
$$

Its purity is

$$
\boxed{
\operatorname{Tr}\rho_e^2
=
\sum_k\lambda_k^2.
}
$$

If one Schmidt coefficient equals one,

$$
\mathcal P=1.
$$

If more than one coefficient is appreciable,

$$
\mathcal P<1.
$$

Thus reduced electronic purity directly detects electron-nuclear entanglement in the
exact closed-system wavefunction.

---

## F. Von Neumann entropy

Diagonalize

$$
\rho_e=V\Lambda V^\dagger
$$

with eigenvalues $\lambda_k$.

Then

$$
\ln\rho_e
=
V(\ln\Lambda)V^\dagger.
$$

Therefore

$$
S_{\rm vN}
=
-\operatorname{Tr}(\rho_e\ln\rho_e)
$$

becomes

$$
\boxed{
S_{\rm vN}
=
-\sum_k
\lambda_k\ln\lambda_k.
}
$$

Terms with $\lambda_k=0$ contribute zero by continuity.

---

## G. Standard error of the mean

For independent samples $X_s$ with variance $\sigma^2$,

$$
\bar X
=
\frac1N\sum_sX_s.
$$

Therefore

$$
\operatorname{Var}(\bar X)
=
\frac{1}{N^2}
\sum_s\operatorname{Var}(X_s)
=
\frac{\sigma^2}{N}.
$$

Hence

$$
\boxed{
\operatorname{SE}(\bar X)
=
\frac{\sigma}{\sqrt N}.
}
$$

The implementation estimates $\sigma$ with the unbiased sample standard deviation.

---

## H. Observed refinement order

Assume

$$
P(h)
=
P^*
+
Ch^p
+
\mathcal O(h^{p+1}).
$$

For refinement ratio $r$,

$$
P(h)-P(h/r)
\approx
Ch^p(1-r^{-p}),
$$

and

$$
P(h/r)-P(h/r^2)
\approx
C(h/r)^p(1-r^{-p}).
$$

Taking the ratio,

$$
\frac{e_h}{e_{h/r}}
\approx
r^p.
$$

Therefore

$$
\boxed{
p
\approx
\frac{\ln(e_h/e_{h/r})}{\ln r}.
}
$$

---

## I. Why sensitivity proxies are not additive

Suppose the method output is

$$
P=P(h,s,b,\ldots),
$$

where $h$ is timestep, $s$ an SPA choice, and $b$ a basis parameter.

Finite differences such as

$$
\Delta_hP
$$

and

$$
\Delta_sP
$$

are evaluated around one shared operating point.

If the parameters interact,

$$
\frac{\partial^2P}{\partial h\partial s}\neq0.
$$

Then the total discrepancy is not

$$
\sqrt{
|\Delta_hP|^2+
|\Delta_sP|^2+\cdots
}.
$$

That formula would require independent random uncertainties.

v0.10 therefore calls these quantities **sensitivity proxies** and reports them
side-by-side.

---

## J. Integrated coupling action and repeated spawning

For parent $i$ and target $a$,

$$
\eta_{ia}(t)
=
|
\dot q_i(t)\cdot d_{ia}(q_i(t))
|.
$$

The accumulated action is

$$
\boxed{
\mathcal A_{ia}(t_1,t_2)
=
\int_{t_1}^{t_2}\eta_{ia}(t)\,dt.
}
$$

Discretely,

$$
\mathcal A
\approx
\sum_n\eta_n\Delta t.
$$

After a spawn, v0.10 may reset

$$
\mathcal A\rightarrow0
$$

and permit later reaccumulation.

A new child is accepted only if the overlap blocker also permits it.

Thus repeated spawning requires both:

$$
\boxed{
\mathcal A\ge\mathcal A_{\rm spawn}
}
$$

and

$$
\boxed{
|S_{\rm child,existing}|<S_{\rm block}.
}
$$
