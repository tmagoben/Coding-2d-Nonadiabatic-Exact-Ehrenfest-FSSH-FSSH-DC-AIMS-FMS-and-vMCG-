# Theoretical Background

## Adiabatic Representation

A molecular Hamiltonian in the absence of an external field is constructed as the sum of electronic and nuclear Hamiltonians, with the addition of an electron-nuclear Coulomb potential

$$
\begin{aligned}
H(r,R)
&= \overbrace{H_e(r)}^{T_e(r)+V_e(r)}
+ \overbrace{H_N(R)}^{T_N(R)+V_N(R)}
+ V_{eN}(r,R), \\
T_e(r)
&= -\sum_i^{n_e}\frac{\nabla_{e,i}^2}{2}, \\
V_e(r)
&= \sum_i^{n_e}\sum_{j>i}^{n_e}\frac{1}{|r_j-r_i|}, \\
T_N(R)
&= -\sum_i^{n_N}\frac{\nabla_{N,i}^2}{2M_i}, \\
V_N(R)
&= \sum_i^{n_N}\sum_{j>i}^{n_N}\frac{Z_jZ_i}{|R_j-R_i|}, \\
V_{eN}(r,R)
&= -\sum_i^{n_N}\sum_j^{n_e}\frac{Z_i}{|r_j-R_i|}.
\end{aligned}
$$

$T_e(r)$, $V_e(r)$, $T_N(R)$, $V_N(R)$, and $V_{eN}(r, R)$ are the electronic kinetic and potential operators, nuclear kinetic and potential operators, and electronic and nuclear interaction potential operator, respectively. The variables $r$ and $R$ denote the electronic and nuclear degrees of freedom, while $n_e$ and $n_N$ are the number of electrons and nuclei in the system. The quantities $M_i$ and $Z_i$ are the nuclear masses and charges.

This molecular Hamiltonian satisfies the time-independent Schr\"odinger equation (TISE),

$$
H(r,R)\Psi_k(r,R)=E_k\Psi_k(r,R).
$$

In the equation above, $\Psi_k(r,R)$ are the eigenfunctions of $H(r,R)$. When $H(r,R)$ is substituted into the time-independent Schr\"odinger equation, the result consists of coupled electronic-nuclear second-order differential equations. At the moment, it is unfeasible to decouple the nuclear and electronic Hilbert spaces by separation of variables. Thus, to decouple the electronic from nuclear degrees of freedom, approximative techniques are required. *(Levine, 2013)*

One technique is to invoke the adiabatic theorem, where the nuclear time scale is considered static relative to the electronic counterpart. Here, the Born--Oppenheimer representation is employed, where one introduces a parameter for the nuclear degrees of freedom in the electronic Hamiltonian, such that *(Born and Oppenheimer)*

$$
H_e(r|R)=T_e(r)+V_e(r)+V_N(R)+V_{eN}(r|R).
$$

The bar in $H_e(r|R)$ represents the parametric dependence of the electronic degrees of freedom on the nuclear coordinates. Because of this parametric dependence, we can obtain eigenstates and eigenfunctions of the electronic Hamiltonian such that

$$
H_e(r|R)\phi_l(r|R)=\epsilon_l(R)\phi_l(r|R).
$$

Here, $\phi_l(r|R)$ are the adiabatic electronic states, and $\epsilon_l(R)$ are the electronic potential energy surfaces. This allows for a redefinition of the molecular Hamiltonian as

$$
H(r,R)=T_N(R)+H_e(r|R).
$$

Next, eigenfunctions of the parameterized molecular Hamiltonian are written using the Born--Huang expansion,

$$
\Psi_k(r,R)=\sum_l \chi_{lk}(R)\phi_l(r|R),
$$

where $\chi_{lk}(R)$ are the nuclear expansion coefficients. *(Born and Huang)* This ansatz is substituted into the redefined molecular Hamiltonian,

$$
\begin{aligned}
H(r,R)\Psi_k(r,R)
&=\sum_l T_N\left[\chi_{lk}(R)\phi_l(r|R)\right]
+H_e(r|R)\chi_{lk}(R)\phi_l(r|R) \\
&=\sum_l T_N\left[\chi_{lk}(R)\phi_l(r|R)\right]
+\epsilon_l(R)\chi_{lk}(R)\phi_l(r|R).
\end{aligned}
$$

To obtain the non-adiabatic elements, we project $\phi_m(r|R)^*$ onto the nuclear kinetic operator. Because $T_N$ contains derivatives with respect to the nuclear coordinates, it acts on both the nuclear coefficient $\chi_{lk}(R)$ and the parametrically nuclear-dependent electronic basis function $\phi_l(r|R)$. Thus,

$$
T_N\left[\chi_{lk}(R)\phi_l(r|R)\right]
=
-\sum_i\frac{1}{2M_i}\nabla_{N,i}^2
\left[\chi_{lk}(R)\phi_l(r|R)\right].
$$

Using the product rule,

$$
\begin{aligned}
\nabla_{N,i}^2
\left[\chi_{lk}(R)\phi_l(r|R)\right]
=
&\left(\nabla_{N,i}^2\chi_{lk}(R)\right)\phi_l(r|R) \\
&+2\left(\nabla_{N,i}\chi_{lk}(R)\right)\cdot
\left(\nabla_{N,i}\phi_l(r|R)\right) \\
&+\chi_{lk}(R)\nabla_{N,i}^2\phi_l(r|R).
\end{aligned}
$$

Projecting onto the electronic state $\phi_m(r|R)$ and inserting a resolution of identity over the electronic adiabatic states gives

$$
\begin{aligned}
\left\langle \phi_m(R) \middle| T_N \middle| \phi_l(R) \right\rangle_r
\chi_{lk}(R)
=
-\sum_i\frac{1}{2M_i}
\bigg[
&\delta_{ml}\nabla_{N,i}^2
+2\tau^{(i)}_{ml}(R)\cdot\nabla_{N,i} \\
&+\nabla_{N,i}\cdot\tau^{(i)}_{ml}(R)
+\sum_n \tau^{(i)}_{mn}(R)\cdot\tau^{(i)}_{nl}(R)
\bigg]\chi_{lk}(R).
\end{aligned}
$$

Here, the anti-Hermitian first-order non-adiabatic coupling vectors are

$$
\tau^{(i)}_{ml}(R)
=
\left\langle \phi_m(R) \middle| \nabla_{N,i}\phi_l(R) \right\rangle_r.
$$

The product term

$$
-\sum_i\frac{1}{2M_i}\sum_n
\tau^{(i)}_{mn}(R)\cdot\tau^{(i)}_{nl}(R)
$$

contributes to the diagonal Born--Oppenheimer correction (DBOC) when $m=l$.

A useful transformation of $\tau^{(i)}_{ml}(R)$ invokes the Hellmann--Feynman theorem,

$$
\nabla_R
\left\langle \phi_m(R) \middle| H_e(R) \middle| \phi_l(R) \right\rangle_r
=
\nabla_R\left[\delta_{ml}\epsilon_m(R)\right],
$$

and substitution of the electronic Hamiltonian $H_e(R)$. This allows one to observe the first-order non-adiabatic coupling vector's reliance on the energy gap between two electronic potential energy surfaces. For $m\neq l$,

$$
\nabla_R
\left\langle \phi_m(R) \middle| H_e(R) \middle| \phi_l(R) \right\rangle_r=0.
$$

Expanding the derivative gives

$$
\begin{aligned}
&\left\langle \nabla_R\phi_m(R) \middle| H_e(R) \middle| \phi_l(R) \right\rangle_r
+\left\langle \phi_m(R) \middle| \nabla_R H_e(R) \middle| \phi_l(R) \right\rangle_r \\
&\qquad
+\left\langle \phi_m(R) \middle| H_e(R) \middle| \nabla_R\phi_l(R) \right\rangle_r=0.
\end{aligned}
$$

Using $H_e(R)\phi_l(R)=\epsilon_l(R)\phi_l(R)$ and $\langle\phi_m(R)|H_e(R)=\epsilon_m(R)\langle\phi_m(R)|$, this becomes

$$
\epsilon_l(R)
\left\langle \nabla_R\phi_m(R) \middle| \phi_l(R) \right\rangle_r
+
\left\langle \phi_m(R) \middle| \nabla_R H_e(R) \middle| \phi_l(R) \right\rangle_r
+
\epsilon_m(R)
\left\langle \phi_m(R) \middle| \nabla_R\phi_l(R) \right\rangle_r
=0.
$$

Rearranging and replacing $\nabla_R$ with $\nabla_{N,i}$, which denotes derivatives with respect to individual nuclear coordinates, gives

$$
\tau^{(i)}_{ml}(R)
=
\left\langle \phi_m(R) \middle| \nabla_{N,i}\phi_l(R) \right\rangle_r
=
\frac{
\left\langle \phi_m(R) \middle| \nabla_{N,i}H_e(R) \middle| \phi_l(R) \right\rangle_r
}{
\epsilon_l(R)-\epsilon_m(R)
}.
$$

*(See discussion of nonadiabatic dynamical simulations near conical intersections via surface-hopping methods.)*

The adiabatic representation accurately describes quantum molecular dynamics in regions absent from potential energy surface intersections. This is a consequence of the presence of non-adiabatic coupling terms $\tau^{(i)}_{mn}$, which become singular in these regions due to their inverse electronic energy gap dependence, $(\epsilon_n(R)-\epsilon_m(R))^{-1}$. Thus, calculating the $\tau^{(i)}_{mn}$ terms requires advanced analytical gradient techniques, which scale steeply with the number of orbitals used in the electronic structure method. *(TNAC; analytical NAC references)* Consequently, one generally searches for a different ansatz to expand the molecular Hamiltonian in.

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

with all non-adiabatic effects condensed into the operator $\Lambda_{ml}(R)$,

$$
\boxed{
\Lambda_{ml}(R)
=
-\sum_i\frac{1}{2M_i}
\left[
2\tau^{(i)}_{ml}(R)\cdot\nabla_{N,i}
+
\nabla_{N,i}\cdot\tau^{(i)}_{ml}(R)
+
\sum_n\tau^{(i)}_{mn}(R)\cdot\tau^{(i)}_{nl}(R)
\right].
}
$$

If all $\tau^{(i)}_{mn}(R)$ are negligible, meaning that electronic potential energy differences are sufficiently large, one may ignore $\Lambda_{ml}(R)$. By doing so, the Born--Oppenheimer equation is realized:

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

As the $\tau_{ml}$ elements are divergent in regions of quasi or exact degeneracy, the question arises whether one can rotate the Hamiltonian to a basis where the derivative couplings $\tau^{(i)}_{mn}=0$. If possible, all non-adiabatic effects would be encapsulated in the potential matrix. *(Beyond Born--Oppenheimer)* As Mead and Truhlar proved in the 1980s, this is generally not possible globally, but the creation of so-called *quasi-diabatic states* is possible, where the derivative couplings $\tau^{(i)}_{mn}$ are minimized. *(Conditions for the definition of a strictly diabatic electronic basis for molecular systems)*

The most trivial method to achieve a diabatic representation from the adiabatic representation is to introduce the *crude-adiabatic basis*,

$$
\Psi_k^{\mathrm{cr}}(r,R)
=
\sum_l\chi_{lk}(R)\phi_l(r|R_0),
$$

where $R_0$ is a reference nuclear geometry. Generally, $R_0$ is chosen to be the minimum of the conical intersection seam, because $R_0$ is a fixed point. For this ansatz, the nuclear kinetic term $T_N$ no longer acts on the electronic eigenfunctions $\phi_l(r|R_0)$.

Next, substituting $\Psi_k^{\mathrm{cr}}(r,R)$ into the molecular TISE gives

$$
\begin{aligned}
H(r,R)\Psi_k^{\mathrm{cr}}(r,R)
&=
\sum_l \phi_l(r|R_0)T_N\chi_{lk}(R)
+
\sum_l \chi_{lk}(R)H_e(r|R)\phi_l(r|R_0) \\
&=E_k\sum_l\chi_{lk}(R)\phi_l(r|R_0).
\end{aligned}
$$

Projecting both sides of the equation onto $\phi_m(r|R_0)^*$, one obtains

$$
\sum_l
\left[
T_N\delta_{ml}
+
\left\langle \phi_m(R_0) \middle| H_e(R) \middle| \phi_l(R_0) \right\rangle_r
\right]
\chi_{lk}(R)
=
E_k\chi_{mk}(R).
$$

Here,

$$
\left\langle \phi_m(R_0) \middle| H_e(R) \middle| \phi_l(R_0) \right\rangle_r
$$

are matrix elements of the electronic Hamiltonian for a fixed reference nuclear geometry. Setting

$$
V(r|R)\equiv V_e(r)+V_{eN}(r|R)+V_N(R),
$$

we define the electronic Hamiltonian as

$$
H_e(r|R)=T_e(r)+V(r|R)
=H_e(r|R_0)-V(r|R_0)+V(r|R).
$$

Rearranging gives a useful representation of the electronic kinetic energy operator,

$$
T_e(r)=H_e(r|R_0)-V(r|R_0).
$$

This allows for the separation of the electronic Hamiltonian $H_e(R_0)$ component from the fixed electronic potential energy $V(r|R_0)$. Evaluating the electronic Hamiltonian matrix elements gives

$$
\begin{aligned}
\left\langle \phi_m(R_0) \middle| H_e(R) \middle| \phi_l(R_0) \right\rangle_r
&=
\left\langle \phi_m(R_0) \middle| H_e(R_0)-V(R_0)+V(R) \middle| \phi_l(R_0) \right\rangle_r \\
&=
\epsilon_l(R_0)\delta_{ml}
+
\left\langle \phi_m(R_0) \middle| V(R)-V(R_0) \middle| \phi_l(R_0) \right\rangle_r.
\end{aligned}
$$

The residual contributions may be redefined as the diabatic potential matrix,

$$
U_{ml}(R)
=
\left\langle \phi_m(R_0) \middle| V(R)-V(R_0) \middle| \phi_l(R_0) \right\rangle_r.
$$

Substituting this into the crude-adiabatic projected equation, the newly defined diabatic Hamiltonian is

$$
\sum_l
\left[
T_N\delta_{ml}
+
\epsilon_l(R_0)\delta_{ml}
+
U_{ml}(R)
\right]
\chi_{lk}(R)
=
E_k\chi_{mk}(R).
$$

By redefining the ansatz in the crude-adiabatic basis, we have effectively transformed all non-adiabatic coupling effects into the off-diagonal diabatic matrix elements. Because the nuclear kinetic operator in the crude-adiabatic ansatz no longer acts on the electronic eigenfunctions, the singularity present in the derivative coupling expression is absent. An additional feature of the diabatic representation is the reduction in the criterion needed for electronic potential energy surface degeneracy. In the adiabatic representation, degeneracies occur when electronic energies are equal and off-diagonal non-adiabatic coupling elements are zero. In the diabatic representation, only the diabatic energies must be degenerate.

The major drawback of using the crude-adiabatic representation is the reliance on $R_0$. Since the nuclear geometry is fixed within the electronic eigenfunctions, an infinite number of individual electronic reference calculations would be needed to map the dynamics of a single problem, which is practically impossible. *(Conical intersections; Introduction to Quantum Mechanics: Time-Dependent references)*
