To obtain the non-adiabatic elements, we project the nuclear kinetic operator onto the electronic adiabatic basis. Starting from the Born--Huang expansion,

$$
\Psi_k(r,R)=\sum_l \chi_{lk}(R)\phi_l(r;R),
$$

the nuclear kinetic operator acts on both the nuclear coefficient and the nuclear-parametric electronic basis function. Thus,

$$
T_N\left[\chi_{lk}(R)\phi_l(r;R)\right]
=
-\sum_i \frac{1}{2M_i}
\nabla_{N,i}^2
\left[\chi_{lk}(R)\phi_l(r;R)\right].
$$

Using the product rule,

$$
\nabla_{N,i}^2
\left[\chi_{lk}(R)\phi_l(r;R)\right]
=
\left(\nabla_{N,i}^2\chi_{lk}(R)\right)\phi_l(r;R)
+
2\left(\nabla_{N,i}\chi_{lk}(R)\right)\cdot
\left(\nabla_{N,i}\phi_l(r;R)\right)
+
\chi_{lk}(R)\nabla_{N,i}^2\phi_l(r;R).
$$

Projecting onto the electronic state \(\phi_m(r;R)\), one obtains

$$
\begin{aligned}
\left\langle \phi_m(R) \middle| T_N \middle| \phi_l(R) \right\rangle_r
\chi_{lk}(R)
=
-\sum_i \frac{1}{2M_i}
\bigg[
&\delta_{ml}\nabla_{N,i}^2
+
2\tau^{(i)}_{ml}(R)\cdot\nabla_{N,i}  \\
&+
\nabla_{N,i}\cdot\tau^{(i)}_{ml}(R)
+
\sum_n
\tau^{(i)}_{mn}(R)\cdot\tau^{(i)}_{nl}(R)
\bigg]\chi_{lk}(R).
\end{aligned}
$$

Here, the first-order non-adiabatic coupling vector is

$$
\tau^{(i)}_{ml}(R)
=
\left\langle
\phi_m(R)
\middle|
\nabla_{N,i}\phi_l(R)
\right\rangle_r .
$$

The term

$$
\sum_n
\tau^{(i)}_{mn}(R)\cdot\tau^{(i)}_{nl}(R)
$$

arises from inserting a resolution of identity over the electronic adiabatic states into the second derivative coupling. For \(m=l\), the corresponding diagonal contribution is associated with the diagonal Born--Oppenheimer correction.

Using this formalism, the molecular Hamiltonian takes the form

$$
\sum_l
\left[
T_N\delta_{ml}
+
\epsilon_m(R)\delta_{ml}
+
\Lambda_{ml}(R)
\right]
\chi_{lk}(R)
=
E_k\chi_{mk}(R),
$$

where all non-adiabatic effects are contained in the operator

$$
\boxed{
\Lambda_{ml}(R)
=
-\sum_i \frac{1}{2M_i}
\left[
2\tau^{(i)}_{ml}(R)\cdot\nabla_{N,i}
+
\nabla_{N,i}\cdot\tau^{(i)}_{ml}(R)
+
\sum_n
\tau^{(i)}_{mn}(R)\cdot\tau^{(i)}_{nl}(R)
\right].
}
$$

If all \(\tau^{(i)}_{mn}(R)\) are negligible, meaning that electronic potential energy differences are sufficiently large, then one may neglect \(\Lambda_{ml}(R)\). In this limit, the Born--Oppenheimer equation is recovered:

$$
\sum_l
\left[
T_N\delta_{ml}
+
\epsilon_l(R)\delta_{ml}
\right]
\chi_{lk}(R)
=
E_k\chi_{mk}(R).
$$
