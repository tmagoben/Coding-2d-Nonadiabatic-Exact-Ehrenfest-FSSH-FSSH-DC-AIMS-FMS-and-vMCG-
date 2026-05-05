# Theoretical Background

## Adiabatic Representation

A molecular Hamiltonian in the absence of an external field is constructed as the sum of electronic and nuclear Hamiltonians with the addition of an electron-nuclear Coulomb potential:

**Equation (2.1)**

$$
\begin{aligned}
H(r,R) &= \overbrace{H_e(r)}^{T_e(r)+V_e(r)} + \overbrace{H_N(R)}^{T_N(R)+V_N(R)} + V_{eN}(r,R), \\
T_e(r) &= -\sum_i^{n_e}\frac{\nabla_{e,i}^2}{2}, \\
V_e(r) &= \sum_i^{n_e}\sum_{j>i}^{n_e}\frac{1}{|r_j-r_i|}, \\
T_N(R) &= -\sum_i^{n_N}\frac{\nabla_{N,i}^2}{2M_i}, \\
V_N(R) &= \sum_i^{n_N}\sum_{j>i}^{n_N}\frac{Z_jZ_i}{|R_j-R_i|}, \\
V_{eN}(r,R) &= -\sum_i^{n_N}\sum_j^{n_e}\frac{Z_i}{|r_j-R_i|}.
\end{aligned}
$$

$T_e(r)$, $V_e(r)$, $T_N(R)$, $V_N(R)$, and $V_{eN}(r,R)$ are the electronic kinetic and potential operators, nuclear kinetic and potential operators, and electronic-nuclear interaction potential operator, respectively. $r$ and $R$ denote the electronic and nuclear degrees of freedom, while $n_e$ and $n_N$ are the number of electrons and nuclei in the system, and $M_i$ and $Z_i$ are the nuclear masses and charges.

This molecular Hamiltonian satisfies the time-independent Schrödinger equation (TISE):

**Equation (2.2)**

$$
H(r,R)\Psi_k(r,R) = E_k\Psi_k(r,R).
$$

In the equation above, $\Psi_k(r,R)$ are the eigenfunctions of $H(r,R)$. When $H(r,R)$ is substituted into the time-independent Schrödinger equation, Equation (2.2) consists of coupled electronic-nuclear second-order differential equations. At the moment, it is unfeasible to decouple the nuclear and electronic Hilbert spaces by separation of variables. Thus, to decouple the electronic from nuclear degrees of freedom, approximate techniques are required. *(cite: `levine2013quantum`)*

One technique is to invoke the adiabatic theorem, where the nuclear time scale is considered static relative to the electronic counterpart. Here, the Born-Oppenheimer representation is employed, where one introduces a parameter for nuclear degrees of freedom in the electronic Hamiltonian such that 

$$
H_e(r\mid R) = T_e(r) + V_e(r) + V_N(R) + V_{eN}(r\mid R).
$$

The bar in $H_e(r\mid R)$ represents the parametric dependence of the electronic degrees of freedom on the nuclear degrees of freedom. Because of this parametric dependence, we can obtain eigenstates and eigenfunctions of the electronic Hamiltonian such that

$$
H_e(r\mid R)\phi_l(r\mid R) = \epsilon_l(R)\phi_l(r\mid R).
$$

Here, $\phi_l(r\mid R)$ are the adiabatic electronic states, and $\epsilon_l(R)$ are the electronic potential energy surfaces. This allows for a redefinition of the molecular Hamiltonian in Equation (2.2) as

$$
H(r,R)=T_N(R)+H_e(r\mid R).
$$

Next, eigenfunctions of the parameterized molecular Hamiltonian are written using the Born-Huang expansion:

$$
\Psi_k(r,R) = \sum_l \chi_{lk}(R)\phi_l(r\mid R),
$$

where $\chi_{lk}(R)$ are the nuclear expansion coefficients. *(cite: `Born_Huang`)* This ansatz is substituted into the redefined molecular Hamiltonian from Equation (2.2):

**Equation (2.5)**

$$
\begin{aligned}
H(r,R)\Psi_k(r,R)
&= \sum_l T_N\left[\chi_{lk}(R)\phi_l(r\mid R)\right]
 + H_e(r\mid R)\chi_{lk}(R)\phi_l(r\mid R) \\
&= \sum_l T_N\left[\chi_{lk}(R)\phi_l(r\mid R)\right]
 + \epsilon_l(R)\chi_{lk}(R)\phi_l(r\mid R).
\end{aligned}
$$

To obtain the non-adiabatic elements, we project $\phi_m^*(r\mid R)$ onto the nuclear kinetic operators in Equation (2.5). Then, inserting a resolution of identity with respect to the electronic adiabatic states gives

$$
\begin{aligned}
\left\langle \phi_m(R) \middle| T_N \middle| \phi_l(R) \right\rangle_r\chi_{lk}
= -\sum_i \frac{1}{2M_i}\Bigg(&\delta_{ml}\nabla_{N,i}^2
+ \tau^{(i)}_{ml}(R)\nabla_{N,i}
+ \nabla_{N,i}\tau^{(i)}_{ml}(R) \\
&+ \sum_n \tau^{(i)}_{mn}(R)\tau^{(i)}_{nl}(R)\Bigg)\chi_{lk}.
\end{aligned}
$$

Here, the anti-Hermitian first-order non-adiabatic coupling vectors are

$$
\tau_{ml}^{(i)}(R)=\left\langle \phi_m(R) \middle| \nabla_{N,i}\phi_l(R) \right\rangle_r,
$$

while the sum of products

$$
-\sum_n\left(\sum_i\frac{1}{2M_i}\tau^{(i)}_{mn}(R)\tau^{(i)}_{nl}(R)\right)
$$

are the diagonal Born-Oppenheimer corrections (DBOC) when $m=l$.

A useful transformation of $\tau_{ml}^{(i)}$ invokes the Hellmann-Feynman theorem,

$$
\nabla_R\left\langle \phi_m(R) \middle| H_e(R) \middle| \phi_l(R) \right\rangle_r
= \nabla_R\delta_{ml}\epsilon_m(R),
$$

and substitution of the electronic Hamiltonian $H_e(R)$. This allows one to observe the first-order non-adiabatic coupling vector's reliance on the energy-gap difference of the two degenerate electronic energy surfaces:

**Equation (2.8)**

$$
\begin{aligned}
\nabla_R\left\langle \phi_m(R) \middle| H_e(R) \middle| \phi_l(R) \right\rangle_r &= 0, \quad m\neq l, \\
\left\langle \nabla_R\phi_m(R) \middle| H_e(R) \middle| \phi_l(R) \right\rangle_r
&+ \left\langle \phi_m(R) \middle| \nabla_R H_e(R) \middle| \phi_l(R) \right\rangle_r \\
&+ \left\langle \phi_m(R) \middle| H_e(R) \middle| \nabla_R\phi_l(R) \right\rangle_r = 0, \\
\epsilon_l(R)\left\langle \nabla_R\phi_m(R) \middle| \phi_l(R) \right\rangle_r
&+ \left\langle \phi_m(R) \middle| \nabla_R H_e(R) \middle| \phi_l(R) \right\rangle_r \\
&+ \epsilon_m(R)\left\langle \phi_m(R) \middle| \nabla_R\phi_l(R) \right\rangle_r = 0.
\end{aligned}
$$

Rearranging Equation (2.8) and replacing the subscript $\nabla_R$ with $\nabla_{N,i}$, denoting individual nuclear coordinates, gives

**Equation (2.9)**

$$
\tau^{(i)}_{ml}(R)
= \left\langle \phi_m(R) \middle| \nabla_{N,i}\phi_l(R) \right\rangle
= \frac{\left\langle \phi_m(R) \middle| \nabla_{N,i}H_e(R) \middle| \phi_l(R) \right\rangle_r}{\epsilon_l(R)-\epsilon_m(R)}.
$$

*(cite: `On_Nonadiabatic_Dynamical_Simulations_near_Conical_Intersections_via_Surface_hopping_Methods`)*

The adiabatic representation accurately describes quantum molecular dynamics in regions absent from potential energy surface intersections. This is a consequence of the presence of non-adiabatic coupling terms $\tau_{mn}^{(i)}$, which become singular in these regions due to their inverse electronic energy-gap dependence, $(\epsilon_n(R)-\epsilon_m(R))^{-1}$. Thus, to calculate the $\tau_{mn}^{(i)}$ terms, advanced analytical gradient techniques are needed, which scale exponentially with the number of orbitals used in the electronic structure method. *(cite: `TNAC`)* Consequently, one generally searches for a different ansatz to expand the molecular Hamiltonian in. *(cite: `Anal_Nac`)*

Using this formalism, the molecular Hamiltonian takes the form

$$
\sum_l \left[T_N\delta_{ml}+\epsilon_m(R)\delta_{ml}+\Lambda_{ml}(R)\right]\chi_{lk}(R)
= E_k\chi_{mk}(R),
$$

with all non-adiabatic effects condensed into the operator $\Lambda_{ml}$:

$$
\Lambda_{ml}
= -\sum_i\frac{1}{2M_i}\left(
\tau^{(i)}_{ml}(R)\nabla_{N,i}
+ \nabla_{N,i}\tau^{(i)}_{ml}(R)
+ \sum_n \tau^{(i)}_{mn}(R)\tau^{(i)}_{nl}(R)
\right).
$$

If all $\tau^{(i)}_{mn}(R)$ are negligible, implying that electronic potential energy differences are significant, one may ignore $\Lambda_{ml}$. By doing so, the Born-Oppenheimer equation is realized:

$$
\sum_l\left(T_N\delta_{ml}+\epsilon_l(R)\delta_{ml}\right)\chi_{lk}(R)=E_k\chi_{mk}(R).
$$

## Diabatic Representation

As the $\tau_{ml}$ elements are divergent in regions of quasi- or exact degeneracy, the question arises whether one can rotate the Hamiltonian to some basis where the derivative couplings $\tau^{(i)}_{mn}=0$. If possible, all non-adiabatic effects would be encapsulated in the potential matrix. *(cite: `Beyond_Born_Oppenheimer`)* As Mead and Truhlar proved in the 1980s, this is generally not the case, but the creation of so-called *quasi-diabatic states* exists, where the minimization of $\tau^{(i)}_{mn}$ can be done. *(cite: `Conditions_for_the_definition_of_a_strictly_diabatic_electronic_basis_for_molecular_systems`)*

The most trivial method to achieve a diabatic representation from an adiabatic representation is to introduce the *crude-adiabatic basis*:

**Equation (2.13)**

$$
\Psi_k^{\mathrm{cr}}(r,R) = \sum_l \chi_{kl}(R)\phi_l(r\mid R_0),
$$

where $R_0$ is a reference point of nuclear geometry. Generally, $R_0$ is chosen to be the minimum of the conical intersection, as $R_0$ is a fixed point. For this ansatz, the nuclear kinetic term $T_N$ no longer acts on the electronic eigenfunctions $\phi_l(r\mid R_0)$.

Next, substituting $\Psi_k^{\mathrm{cr}}(r,R)$ into Equation (2.2), the molecular Hamiltonian takes the form

**Equation (2.12)**

$$
\begin{aligned}
H(r,R)\Psi_k^{\mathrm{cr}}(r,R)
&= \sum_l \phi_l(r\mid R_0)T_N\chi_{lk}(R)
+ \sum_l \chi_{lk}(R)H_e(r\mid R)\phi_l(r\mid R_0) \\
&= E_k\sum_l\chi_{lk}(R)\phi_l(r\mid R_0).
\end{aligned}
$$

Projecting both sides of the equation onto $\phi_m^*(r\mid R_0)$, one arrives at

**Equation (2.15)**

$$
\sum_l\left[T_N\delta_{ml}
+ \left\langle \phi_m(R_0) \middle| H_e(R) \middle| \phi_l(R_0) \right\rangle_r\right]\chi_{lk}(R)
= E_k\chi_{mk}(R).
$$

Here,

$$
\left\langle \phi_m(R_0) \middle| H_e(R) \middle| \phi_l(R_0) \right\rangle_r
$$

are matrix elements of the electronic Hamiltonian for a fixed nuclear geometry. Setting

$$
V(r\mid R) \equiv V_e(r)+V_{eN}(r\mid R)+V_N(R),
$$

we define the electronic Hamiltonian as

**Equation (diabatic)**

$$
H_e(r\mid R)=T_e(r)+V(r\mid R)=H_e(r\mid R_0)-V(r\mid R_0)+V(r\mid R).
$$

Rearranging the diabatic equation to the form below gives a useful representation of the electronic kinetic energy operator. This allows for the separation of the electronic Hamiltonian $H_e(R_0)$ component from the fixed electronic potential energy $V(r\mid R_0)$:

**Equation (2.16)**

$$
\begin{aligned}
T_e(r) &= H_e(r\mid R_0)-V(r\mid R_0), \\
&= T_e(r)+V_e(r)+V_{eN}(r\mid R_0)+V_N(R)
-\left[V_e(r)+V_{eN}(r\mid R_0)+V_N(R)\right].
\end{aligned}
$$

Evaluating the electronic Hamiltonian matrix elements, the fixed electronic potential energy surfaces are obtained in addition to a matrix containing the residual potential energy contributions:

$$
\begin{aligned}
\left\langle \phi_m(R_0) \middle| H_e(R) \middle| \phi_l(R_0) \right\rangle_r
&= \left\langle \phi_m(R_0) \middle| H_e(R_0)-V(R_0)+V(R) \middle| \phi_l(R_0) \right\rangle_r \\
&= \epsilon_l(R_0)\delta_{ml}
+ \left\langle \phi_m(R_0) \middle| V(R)-V(R_0) \middle| \phi_l(R_0) \right\rangle_r.
\end{aligned}
$$

The residual contributions may be redefined as the diabatic potential matrix:

**Equation (2.18)**

$$
U_{ml}=\left\langle \phi_m(R_0) \middle| V(R)-V(R_0) \middle| \phi_l(R_0) \right\rangle_r.
$$

Substituting Equation (2.18) into Equation (2.15), the newly defined diabatic Hamiltonian is *(cite: `Introduction_to_Quantum_Mechanics-Time-Dependent`)*

**Equation (2.19)**

$$
\sum_l\left[T_N\delta_{ml}+\epsilon_m(R_0)\delta_{ml}+U_{ml}\right]\chi_{lk}(R)
= E_k\chi_{mk}(R).
$$

By redefining our ansatz using the crude-adiabatic basis, we have effectively transformed all non-adiabatic coupling effects into the off-diagonal diabatic matrix elements. Because the nuclear kinetic operator in Equation (2.13) no longer acts on the electronic eigenfunctions, the singularity present in Equation (2.9) is absent. An additional feature of the diabatic representation is the reduction in the criterion needed for electronic potential energy surface degeneracy. In the adiabatic representation, degeneracies occur when electronic energies are equal and off-diagonal non-adiabatic coupling elements are zero. In the diabatic representation, only the energies must be degenerate.

The major drawback of using the crude-adiabatic representation is the reliance on $R_0$. Because we have fixed our nuclear geometry within the electronic eigenfunctions, an infinite number of individual electronic reference calculations would be needed to map the dynamics of a single problem, which is practically impossible. *(cite: `Conical_Intersections`; `Introduction_to_Quantum_Mechanics-Time-Dependent`)*
