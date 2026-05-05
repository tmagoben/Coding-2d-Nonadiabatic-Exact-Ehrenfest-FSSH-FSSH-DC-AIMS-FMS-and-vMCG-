# Theoretical Background

## Adiabatic Representation

A molecular Hamiltonian in the absence of an external field is constructed as the sum of electronic and nuclear Hamiltonians, with the addition of the electron--nuclear Coulomb potential:

$$
\begin{aligned}
H(r,R) &= H_e(r)+H_N(R)+V_{eN}(r,R), \\
H_e(r) &= T_e(r)+V_e(r), \\
H_N(R) &= T_N(R)+V_N(R), \\
T_e(r) &= -\sum_{i=1}^{n_e}\frac{\nabla_{e,i}^2}{2}, \\
V_e(r) &= \sum_{i=1}^{n_e}\sum_{j>i}^{n_e}\frac{1}{\lvert r_j-r_i\rvert}, \\
T_N(R) &= -\sum_{i=1}^{n_N}\frac{\nabla_{N,i}^2}{2M_i}, \\
V_N(R) &= \sum_{i=1}^{n_N}\sum_{j>i}^{n_N}\frac{Z_jZ_i}{\lvert R_j-R_i\rvert}, \\
V_{eN}(r,R) &= -\sum_{i=1}^{n_N}\sum_{j=1}^{n_e}\frac{Z_i}{\lvert r_j-R_i\rvert}.
\end{aligned}
$$

Here, $T_e(r)$, $V_e(r)$, $T_N(R)$, $V_N(R)$, and $V_{eN}(r,R)$ are the electronic kinetic operator, electronic potential operator, nuclear kinetic operator, nuclear potential operator, and electron--nuclear interaction potential operator, respectively. The variables $r$ and $R$ denote the electronic and nuclear degrees of freedom, while $n_e$ and $n_N$ are the numbers of electrons and nuclei. The quantities $M_i$ and $Z_i$ are the mass and charge of nucleus $i$.

This molecular Hamiltonian satisfies the time-independent Schrödinger equation:

$$
H(r,R)\Psi_k(r,R)=E_k\Psi_k(r,R).
$$

In the equation above, $\Psi_k(r,R)$ are the eigenfunctions of $H(r,R)$. When $H(r,R)$ is substituted into the time-independent Schrödinger equation, one obtains coupled electronic--nuclear second-order differential equations. In general, it is not feasible to decouple the nuclear and electronic Hilbert spaces by separation of variables. Thus, approximate techniques are required.

One technique is to invoke the adiabatic theorem, where the nuclear time scale is considered slow relative to the electronic time scale. In the Born--Oppenheimer representation, the nuclear degrees of freedom enter the electronic Hamiltonian as parameters:

$$
H_e(r;R)=T_e(r)+V_e(r)+V_N(R)+V_{eN}(r;R).
$$

The semicolon in $H_e(r;R)$ denotes the parametric dependence of the electronic Hamiltonian on the nuclear coordinates. Because of this parametric dependence, one obtains eigenstates and eigenvalues of the electronic Hamiltonian:

$$
H_e(r;R)\phi_l(r;R)=\epsilon_l(R)\phi_l(r;R).
$$

Here, $\phi_l(r;R)$ are the adiabatic electronic states and $\epsilon_l(R)$ are the electronic potential energy surfaces. The molecular Hamiltonian may then be written as

$$
H(r,R)=T_N(R)+H_e(r;R).
$$

The eigenfunctions of the parametrized molecular Hamiltonian are written using the Born--Huang expansion:

$$
\Psi_k(r,R)=\sum_l \chi_{lk}(R)\phi_l(r;R),
$$

where $\chi_{lk}(R)$ are nuclear expansion coefficients. Substituting this ansatz into the molecular Schrödinger equation gives

$$
\begin{aligned}
H(r,R)\Psi_k(r,R)
&=\sum_l T_N\big[\chi_{lk}(R)\phi_l(r;R)\big]
+\sum_l H_e(r;R)\chi_{lk}(R)\phi_l(r;R) \\
&=\sum_l T_N\big[\chi_{lk}(R)\phi_l(r;R)\big]
+\sum_l \epsilon_l(R)\chi_{lk}(R)\phi_l(r;R).
\end{aligned}
$$

To obtain the non-adiabatic coupling terms, project the nuclear kinetic operator onto the electronic adiabatic basis. Since $T_N$ differentiates with respect to nuclear coordinates, it acts on both the nuclear coefficient $\chi_{lk}(R)$ and the nuclear-parametric electronic basis function $\phi_l(r;R)$:

$$
T_N\!\left[\chi_{lk}(R)\phi_l(r;R)\right]
=
-\sum_i \frac{1}{2M_i}
\nabla_{N,i}^{2}
\!\left[\chi_{lk}(R)\phi_l(r;R)\right].
$$

Using the product rule,

$$
\begin{aligned}
\nabla_{N,i}^2\left[\chi_{lk}(R)\phi_l(r;R)\right]
&=\left(\nabla_{N,i}^2\chi_{lk}(R)\right)\phi_l(r;R) \\
&\quad +2\left(\nabla_{N,i}\chi_{lk}(R)\right)\cdot\left(\nabla_{N,i}\phi_l(r;R)\right) \\
&\quad +\chi_{lk}(R)\nabla_{N,i}^2\phi_l(r;R).
\end{aligned}
$$

Projecting onto $\phi_m(r;R)$ gives the nuclear kinetic energy matrix element:

$$
\begin{aligned}
\left\langle \phi_m(R) \right\rvert T_N \left\lvert \phi_l(R) \right\rangle_r\chi_{lk}(R)
= -\sum_i\frac{1}{2M_i}
\bigg[
&\delta_{ml}\nabla_{N,i}^2
+2\tau^{(i)}_{ml}(R)\cdot\nabla_{N,i} \\
&+\nabla_{N,i}\cdot\tau^{(i)}_{ml}(R)
+\sum_n \tau^{(i)}_{mn}(R)\cdot\tau^{(i)}_{nl}(R)
\bigg]\chi_{lk}(R).
\end{aligned}
$$

The first-order non-adiabatic coupling vector is

$$
\tau^{(i)}_{ml}(R)
=
\left\langle \phi_m(R) \middle\vert \nabla_{N,i}\phi_l(R) \right\rangle_r.
$$

The product term

$$
\sum_n \tau^{(i)}_{mn}(R)\cdot\tau^{(i)}_{nl}(R)
$$

arises from inserting a resolution of identity over the electronic adiabatic states into the second derivative coupling. For $m=l$, the corresponding diagonal contribution is associated with the diagonal Born--Oppenheimer correction.

A useful transformation of $\tau^{(i)}_{ml}$ follows from the Hellmann--Feynman relation. For $m\neq l$,

$$
\nabla_R\left\langle\phi_m(R)\right\rvert H_e(R)\left\lvert\phi_l(R)\right\rangle_r=0.
$$

Expanding this derivative gives

$$
\begin{aligned}
&\left\langle\nabla_R\phi_m(R)\right\rvert H_e(R)\left\lvert\phi_l(R)\right\rangle_r
+\left\langle\phi_m(R)\right\rvert\nabla_R H_e(R)\left\lvert\phi_l(R)\right\rangle_r \\
&\quad +\left\langle\phi_m(R)\right\rvert H_e(R)\left\lvert\nabla_R\phi_l(R)\right\rangle_r=0.
\end{aligned}
$$

Using the electronic eigenvalue equation,

$$
\epsilon_l(R)\left\langle\nabla_R\phi_m(R)\middle\vert\phi_l(R)\right\rangle_r
+\left\langle\phi_m(R)\right\rvert\nabla_R H_e(R)\left\lvert\phi_l(R)\right\rangle_r
+\epsilon_m(R)\left\langle\phi_m(R)\middle\vert\nabla_R\phi_l(R)\right\rangle_r=0.
$$

For $m\neq l$, orthonormality implies

$$
\left\langle\nabla_R\phi_m(R)\middle\vert\phi_l(R)\right\rangle_r
=
-\left\langle\phi_m(R)\middle\vert\nabla_R\phi_l(R)\right\rangle_r.
$$

Therefore,

$$
\tau^{(i)}_{ml}(R)
=
\left\langle\phi_m(R)\middle\vert\nabla_{N,i}\phi_l(R)\right\rangle_r
=
\frac{
\left\langle\phi_m(R)\right\rvert\nabla_{N,i}H_e(R)\left\lvert\phi_l(R)\right\rangle_r
}{
\epsilon_l(R)-\epsilon_m(R)
}.
$$

The adiabatic representation accurately describes quantum molecular dynamics in regions away from potential energy surface intersections. This is because the non-adiabatic coupling vectors $\tau^{(i)}_{mn}(R)$ become large near degeneracies through their inverse dependence on the electronic energy gap $\epsilon_n(R)-\epsilon_m(R)$.

Using this formalism, the molecular Hamiltonian takes the form

$$
\sum_l
\left[
T_N\delta_{ml}
+\epsilon_m(R)\delta_{ml}
+\Lambda_{ml}(R)
\right]
\chi_{lk}(R)
=
E_k\chi_{mk}(R),
$$

where all non-adiabatic effects are contained in the operator

$$
\Lambda_{ml}(R)
=
-\sum_i\frac{1}{2M_i}
\left[
2\tau^{(i)}_{ml}(R)\cdot\nabla_{N,i}
+\nabla_{N,i}\cdot\tau^{(i)}_{ml}(R)
+\sum_n\tau^{(i)}_{mn}(R)\cdot\tau^{(i)}_{nl}(R)
\right].
$$

If all $\tau^{(i)}_{mn}(R)$ are negligible, meaning that electronic potential energy differences are sufficiently large, one may ignore $\Lambda_{ml}(R)$. By doing so, the Born--Oppenheimer equation is recovered:

$$
\sum_l
\left[
T_N\delta_{ml}+\epsilon_l(R)\delta_{ml}
\right]
\chi_{lk}(R)
=
E_k\chi_{mk}(R).
$$

## Diabatic Representation

As the $\tau_{ml}$ elements are divergent in regions of quasi-degeneracy or exact degeneracy, the question arises whether one can rotate the Hamiltonian to a basis where the derivative couplings satisfy

$$
\tau^{(i)}_{mn}=0.
$$

If this were possible, all non-adiabatic effects would be encapsulated in the potential matrix. In general, a strictly diabatic basis does not exist globally, but quasi-diabatic states can be constructed by minimizing the derivative couplings.

The most direct way to construct a diabatic representation from an adiabatic one is to introduce the crude adiabatic basis:

$$
\Psi_k^{\mathrm{cr}}(r,R)=\sum_l\chi_{lk}(R)\phi_l(r;R_0),
$$

where $R_0$ is a fixed reference nuclear geometry. Often, $R_0$ is chosen near the minimum of the conical intersection seam. For this ansatz, the nuclear kinetic operator $T_N$ no longer acts on the electronic basis functions $\phi_l(r;R_0)$.

Substituting $\Psi_k^{\mathrm{cr}}(r,R)$ into the molecular Schrödinger equation gives

$$
\begin{aligned}
H(r,R)\Psi_k^{\mathrm{cr}}(r,R)
&=
\sum_l\phi_l(r;R_0)T_N\chi_{lk}(R)
+
\sum_l\chi_{lk}(R)H_e(r;R)\phi_l(r;R_0) \\
&=E_k\sum_l\chi_{lk}(R)\phi_l(r;R_0).
\end{aligned}
$$

Projecting both sides onto $\phi_m(r;R_0)$ gives

$$
\sum_l
\left[
T_N\delta_{ml}
+
\left\langle\phi_m(R_0)\right\rvert H_e(R)\left\lvert\phi_l(R_0)\right\rangle_r
\right]
\chi_{lk}(R)
=
E_k\chi_{mk}(R).
$$

Here,

$$
\left\langle\phi_m(R_0)\right\rvert H_e(R)\left\lvert\phi_l(R_0)\right\rangle_r
$$

are matrix elements of the electronic Hamiltonian in a fixed electronic basis. Defining

$$
V(r;R)=V_e(r)+V_{eN}(r;R)+V_N(R),
$$

the electronic Hamiltonian may be written as

$$
H_e(r;R)=T_e(r)+V(r;R)=H_e(r;R_0)-V(r;R_0)+V(r;R).
$$

Therefore,

$$
T_e(r)=H_e(r;R_0)-V(r;R_0).
$$

Evaluating the electronic Hamiltonian matrix elements in the crude adiabatic basis gives

$$
\begin{aligned}
\left\langle\phi_m(R_0)\right\rvert H_e(R)\left\lvert\phi_l(R_0)\right\rangle_r
&=
\left\langle\phi_m(R_0)\right\rvert
H_e(R_0)-V(R_0)+V(R)
\left\lvert\phi_l(R_0)\right\rangle_r \\
&=\epsilon_l(R_0)\delta_{ml}
+
\left\langle\phi_m(R_0)\right\rvert
V(R)-V(R_0)
\left\lvert\phi_l(R_0)\right\rangle_r.
\end{aligned}
$$

The residual contribution defines the diabatic potential matrix:

$$
U_{ml}(R)
=
\left\langle\phi_m(R_0)\right\rvert
V(R)-V(R_0)
\left\lvert\phi_l(R_0)\right\rangle_r.
$$

Substituting this definition into the crude-adiabatic projected equation gives the diabatic Hamiltonian:

$$
\sum_l
\left[
T_N\delta_{ml}
+
\epsilon_m(R_0)\delta_{ml}
+
U_{ml}(R)
\right]
\chi_{lk}(R)
=
E_k\chi_{mk}(R).
$$

By redefining the ansatz using the crude adiabatic basis, the non-adiabatic coupling effects are moved into the off-diagonal diabatic potential matrix elements. Because the nuclear kinetic operator no longer acts on the electronic basis functions, the singular derivative coupling term present in the adiabatic representation is absent.

An additional feature of the diabatic representation is the reduction in the criterion needed for electronic potential energy surface degeneracy. In the adiabatic representation, degeneracies require equality of electronic energies and appropriate coupling conditions. In the diabatic representation, degeneracy is expressed directly through the diabatic potential matrix.

The major drawback of using the crude adiabatic representation is its reliance on the fixed reference geometry $R_0$. Because the electronic basis functions are fixed at a single nuclear geometry, many reference calculations would be required to map the dynamics of a general molecular problem, which is usually impractical.
