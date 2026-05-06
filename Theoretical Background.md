# Theoretical Background

This section briefly surveys the adiabatic and diabatic representations, emphasizing their role in simulations of nonadiabatic molecular dynamics. A general background in quantum mechanics is assumed. Unless otherwise stated, atomic units are used, so that \(\hbar=m_e=e=1\).

## Adiabatic Representation

For a molecule in the absence of an external field, the molecular Hamiltonian can be written as the sum of electronic kinetic energy, nuclear kinetic energy, electron-electron repulsion, nucleus-nucleus repulsion, and electron-nucleus attraction terms:

$$
\begin{array}{rl}
H(r,R) &= T_e(r)+T_N(R)+V_{ee}(r)+V_{NN}(R)+V_{eN}(r,R), \
T_e(r) &= -\sum_{i=1}^{n_e}\frac{1}{2}\nabla_{e,i}^{2}, \
T_N(R) &= -\sum_{A=1}^{n_N}\frac{1}{2M_A}\nabla_{N,A}^{2}, \
V_{ee}(r) &= \sum_{i=1}^{n_e}\sum_{j>i}^{n_e}\frac{1}{|r_j-r_i|}, \
V_{NN}(R) &= \sum_{A=1}^{n_N}\sum_{B>A}^{n_N}\frac{Z_AZ_B}{|R_B-R_A|}, \
V_{eN}(r,R) &= -\sum_{A=1}^{n_N}\sum_{i=1}^{n_e}\frac{Z_A}{|r_i-R_A|}.
\end{array}
$$

Here, $\(r\)$ denotes all electronic coordinates and $\(R\)$ denotes all nuclear coordinates. The quantities $\(n_e\)$ and $\(n_N\)$ are the numbers of electrons and nuclei, while $\(M_A\)$ and $\(Z_A\)$ are the mass and charge of nucleus $\(A\)$. The full molecular eigenvalue problem is

$$
H(r,R)\Psi_k(r,R)=E_k\Psi_k(r,R).
$$

The direct solution of this equation is generally difficult because the electronic and nuclear coordinates are coupled. The Born-Oppenheimer approximation exploits the separation between fast electronic motion and slower nuclear motion by treating the nuclear coordinates as fixed parameters in the electronic eigenvalue problem.

For fixed $\(R\)$, define the electronic Hamiltonian

$$
H_{\mathrm{el}}(r;R)=T_e(r)+V_{ee}(r)+V_{eN}(r,R).
$$

The clamped-nuclei electronic eigenvalue equation is

$$
H_{\mathrm{el}}(r;R)\phi_l(r;R)=E_l^{\mathrm{el}}(R)\phi_l(r;R).
$$

It is often convenient to absorb the nuclear repulsion into the electronic energy surface:

$$
\epsilon_l(R)=E_l^{\mathrm{el}}(R)+V_{NN}(R).
$$

With this convention, the molecular Hamiltonian can be written as

$$
H(r,R)=T_N(R)+H_e(r;R),
$$

where

$$
H_e(r;R)=H_{\mathrm{el}}(r;R)+V_{NN}(R).
$$

The adiabatic electronic states $\(\phi_l(r;R)\)$ form a nuclear-coordinate-dependent electronic basis. The molecular wavefunction may therefore be expanded in the Born-Huang form

$$
\Psi_k(r,R)=\sum_l \chi_{lk}(R)\phi_l(r;R),
$$

where $\(\chi_{lk}(R)\)$ is the nuclear wavefunction amplitude on electronic state $\(l\)$ for molecular eigenstate $\(k\)$.

Substituting the Born-Huang expansion into the molecular Schrödinger equation gives

$$
\begin{array}{rl}
H(r,R)\Psi_k(r,R)
=&
\displaystyle\sum_l T_N[\chi_{lk}(R)\phi_l(r;R)]
+
\sum_l H_e(r;R)\chi_{lk}(R)\phi_l(r;R)
\
=&
\displaystyle\sum_l T_N[\chi_{lk}(R)\phi_l(r;R)]
+
\sum_l \epsilon_l(R)\chi_{lk}(R)\phi_l(r;R).
\end{array}
$$

The nonadiabatic coupling terms arise because $\(T_N\)$ differentiates with respect to nuclear coordinates. Therefore, it acts not only on the nuclear coefficient $\(\chi_{lk}(R)\)$, but also on the nuclear-parametric electronic basis function $\(\phi_l(r;R)\)$:

$$
T_N[\chi_{lk}(R)\phi_l(r;R)]=
-\sum_A \frac{1}{2M_A}
\nabla_{N,A}^{2}[\chi_{lk}(R)\phi_l(r;R)].
$$

Using the product rule,

$$
\begin{array}{rl}
\nabla_{N,A}^{2}[\chi_{lk}(R)\phi_l(r;R)]
=&
(\nabla_{N,A}^{2}\chi_{lk}(R))\phi_l(r;R)
\
&
+2(\nabla_{N,A}\chi_{lk}(R))\cdot(\nabla_{N,A}\phi_l(r;R))
\
&
+\chi_{lk}(R)\nabla_{N,A}^{2}\phi_l(r;R).
\end{array}
$$

Define the first-order derivative coupling vector

$$
\tau^{(A)}_{ml}(R) =\langle \phi_m(r;R)|\nabla_{N,A}\phi_l(r;R)\rangle_r,
$$

where $\(\langle\cdots\rangle_r\)$ denotes integration over electronic coordinates only. Also, define the second-derivative coupling

$$
d^{(A)}_{ml}(R) =
\langle \phi_m(r;R)|\nabla_{N,A}^{2}\phi_l(r;R)\rangle_r.
$$

Projecting the nuclear kinetic operator onto the electronic state $\(\phi_m\)$ gives

$$
\begin{array}{rl}
\langle \phi_m|T_N|\phi_l\rangle_r\chi_{lk}(R)=&\displaystyle-\sum_A\frac{1}{2M_A}\left[\delta_{ml}\nabla_{N,A}^{2}+2\tau^{(A)}_{ml}(R)\cdot\nabla_{N,A}+d^{(A)}_{ml}(R)\right]\chi_{lk}(R).\end{array}
$$

The second-derivative coupling can be rewritten by inserting a resolution of identity over the adiabatic electronic states:

$$
d^{(A)}_{ml}(R) =
\nabla_{N,A}\cdot\tau^{(A)}_{ml}(R)
+
\sum_n
\tau^{(A)}_{mn}(R)\cdot\tau^{(A)}_{nl}(R).
$$

Therefore,

$$
\begin{array}{rl}
\langle \phi_m|T_N|\phi_l\rangle_r\chi_{lk}(R)
=&
\displaystyle
-\sum_A\frac{1}{2M_A}
\left[
\delta_{ml}\nabla_{N,A}^{2}
+
2\tau^{(A)}_{ml}(R)\cdot\nabla_{N,A}
\right.
\
&
\displaystyle
\left.
+
\nabla_{N,A}\cdot\tau^{(A)}_{ml}(R)
+
\sum_n
\tau^{(A)}_{mn}(R)\cdot\tau^{(A)}_{nl}(R)
\right]\chi_{lk}(R).
\end{array}
$$

For $\(m=l\)$, the diagonal part of the second-derivative coupling contributes to the diagonal Born-Oppenheimer correction. In a real adiabatic gauge with $\(\tau^{(A)}_{mm}=0\)$, this correction is commonly written as

$$
E_{\mathrm{DBOC},m}(R) =
\sum_A\frac{1}{2M_A}
\sum_{n\ne m}
|\tau^{(A)}_{mn}(R)|^2.
$$

A useful relation for the off-diagonal derivative coupling follows from differentiating the electronic eigenvalue equation. For $\(m\ne l\)$,

$$
\langle \phi_m(R)|H_e(R)|\phi_l(R)\rangle_r=0.
$$

Differentiating with respect to a nuclear coordinate gives

$$
\begin{array}{l}
\langle \nabla_{N,A}\phi_m(R)|H_e(R)|\phi_l(R)\rangle_r +\langle \phi_m(R)|\nabla_{N,A}H_e(R)|\phi_l(R)\rangle_r\
\qquad
+
\langle \phi_m(R)|H_e(R)|\nabla_{N,A}\phi_l(R)\rangle_r=
0.
\end{array}
$$

Using the electronic eigenvalue equation and the orthonormality relation

$$
\langle \nabla_{N,A}\phi_m(R)|\phi_l(R)\rangle_r=-\langle \phi_m(R)|\nabla_{N,A}\phi_l(R)\rangle_r,\qquad m\ne l,
$$

one obtains

$$
\tau^{(A)}_{ml}(R) =\langle \phi_m(R)|\nabla_{N,A}\phi_l(R)\rangle_r =\frac{\langle \phi_m(R)|\nabla_{N,A}H_e(R)|\phi_l(R)\rangle_r}{\epsilon_l(R)-\epsilon_m(R)},\qquad m\ne l.
$$

This expression is valid away from exact degeneracies and assumes a differentiable choice of adiabatic electronic phases. As $\(\epsilon_l(R)-\epsilon_m(R)\)$ becomes small, the expression becomes ill-conditioned and derivative couplings may become large. Thus, near avoided crossings, conical intersections, or other near-degeneracies, nonadiabatic couplings generally cannot be neglected.

The coupled nuclear equations in the adiabatic representation can be written as

$$
\sum_l\left[T_N\delta_{ml}+\epsilon_l(R)\delta_{ml}+\Lambda_{ml}(R)\right]\chi_{lk}(R)=E_k\chi_{mk}(R),
$$

where all nonadiabatic derivative-coupling effects are collected in

$$
\Lambda_{ml}(R)=-\sum_A\frac{1}{2M_A}\left[2\tau^{(A)}_{ml}(R)\cdot\nabla_{N,A}+\nabla_{N,A}\cdot\tau^{(A)}_{ml}(R)+\sum_n\tau^{(A)}_{mn}(R)\cdot\tau^{(A)}_{nl}(R)\right].
$$

If the derivative couplings $\(\tau^{(A)}_{mn}(R)\)$ are negligible over the dynamically relevant region of nuclear configuration space, then 

$\(\Lambda_{ml}(R)\)$ may be neglected. The uncoupled Born-Oppenheimer nuclear equation is then recovered:

$$
\left[T_N+\epsilon_m(R)\right]\chi_{mk}(R)=E_k\chi_{mk}(R).
$$

Equivalently,

$$
\sum_l\left[T_N\delta_{ml}+\epsilon_l(R)\delta_{ml}\right]\chi_{lk}(R)=E_k\chi_{mk}(R).
$$

This Born-Oppenheimer form underlies much of practical molecular electronic-structure theory. At fixed nuclear geometries, electronic-structure calculations provide potential energy surfaces, nuclear gradients, spin-orbit coupling matrix elements, and, when required, nonadiabatic derivative couplings.

## Diabatic Representation

The adiabatic representation is natural for electronic-structure calculations because it diagonalizes the electronic Hamiltonian at each nuclear geometry. However, it can become inconvenient near degeneracies because derivative couplings may become large or singular. This motivates the search for representations in which nonadiabatic effects are shifted from derivative couplings into off-diagonal potential couplings.

A diabatic representation is a basis in which the first-order derivative couplings vanish or are made small:

$$
\tau^{(A)}_{mn}(R)=\langle \varphi_m(r;R)|\nabla_{N,A}\varphi_n(r;R)\rangle_r\approx 0.
$$

If the derivative couplings were exactly zero, the nuclear kinetic operator would remain diagonal in the electronic basis, and nonadiabatic effects would appear through the off-diagonal elements of the electronic potential matrix. In general, a globally valid strictly diabatic basis cannot be constructed for arbitrary molecular systems because derivative couplings cannot usually be transformed away everywhere in nuclear-coordinate space. Nevertheless, local diabatic, quasi-diabatic, or model diabatic bases are often useful.

A simple illustrative route toward a diabatic-like representation is the crude adiabatic basis. In this construction, the electronic basis functions are fixed at a reference geometry $\(R_0\)$:

$$
\Psi_k^{\mathrm{cr}}(r,R)=\sum_l\chi_{lk}(R)\phi_l(r;R_0).
$$

Because $\(\phi_l(r;R_0)\)$ does not depend on the current nuclear coordinate \(R\), the nuclear kinetic operator acts only on the nuclear coefficients:

$$
T_N[\chi_{lk}(R)\phi_l(r;R_0)] = \phi_l(r;R_0)T_N\chi_{lk}(R).
$$

Substituting the crude adiabatic expansion into the molecular Schrödinger equation gives

$$
\begin{array}{rl}
H(r,R)\Psi_k^{\mathrm{cr}}(r,R)
=&
\displaystyle
\sum_l
\phi_l(r;R_0)T_N\chi_{lk}(R)
+
\sum_l
\chi_{lk}(R)H_e(r;R)\phi_l(r;R_0)
\
=&
\displaystyle
E_k\sum_l
\chi_{lk}(R)\phi_l(r;R_0).
\end{array}
$$

Projecting onto \(\phi_m(r;R_0)\) gives

$$
\sum_l\left[T_N\delta_{ml}+V^{\mathrm{cr}}_{ml}(R)\right]\chi_{lk}(R) = E_k\chi_{mk}(R),
$$

where

$$
V^{\mathrm{cr}}_{ml}(R) = \langle \phi_m(r;R_0)|H_e(r;R)|\phi_l(r;R_0)\rangle_r.
$$

Define

$$
V(r;R)=V_{ee}(r)+V_{eN}(r,R)+V_{NN}(R).
$$

Then

$$
H_e(r;R)=T_e(r)+V(r;R),
$$

and

$$
H_e(r;R)=H_e(r;R_0)-V(r;R_0)+V(r;R).
$$

Therefore, the crude-basis potential matrix becomes

$$
\begin{array}{rl}
V^{\mathrm{cr}}_{ml}(R)
=&
\langle \phi_m(r;R_0)|
H_e(r;R_0)-V(r;R_0)+V(r;R)
|\phi_l(r;R_0)\rangle_r
\
=&
\epsilon_l(R_0)\delta_{ml}
+
\langle \phi_m(r;R_0)|
V(r;R)-V(r;R_0)
|\phi_l(r;R_0)\rangle_r.
\end{array}
$$

The residual contribution defines the diabatic potential matrix

$$
U_{ml}(R)=\langle \phi_m(r;R_0)|V(r;R)-V(r;R_0)|\phi_l(r;R_0)\rangle_r.
$$

The crude-diabatic nuclear equation is therefore

$$
\sum_l\left[T_N\delta_{ml}+\epsilon_l(R_0)\delta_{ml}+U_{ml}(R)\right]\chi_{lk}(R)=E_k\chi_{mk}(R).
$$

This construction illustrates how derivative-coupling effects can be shifted from the nuclear kinetic operator into a potential-energy matrix. In this representation, off-diagonal elements of $\(U_{ml}(R)\)$ couple the nuclear amplitudes associated with different electronic basis functions.

It is important not to confuse diabatic crossings with adiabatic degeneracies. For a two-state diabatic potential matrix

$$
V_{\mathrm{dia}}(R) = \begin{pmatrix} V_{11}(R) & V_{12}(R) \\ 
V_{21}(R) & V_{22}(R) \\ 
\end{pmatrix}
$$

The adiabatic energies are the eigenvalues of this matrix. For a real symmetric two-state model,

$$
E_{\pm}(R) = \frac{V_{11}(R)+V_{22}(R)}{2}\pm\sqrt{\left[\frac{V_{11}(R)-V_{22}(R)}{2}\right]^2+V_{12}(R)^2}.
$$

Thus, equality of the diagonal diabatic potentials, $\(V_{11}(R)=V_{22}(R)\)$, does not by itself imply an adiabatic degeneracy if the off-diagonal coupling $\(V_{12}(R)\)$ is nonzero. Instead, the off-diagonal diabatic coupling controls the size of the avoided crossing after diagonalization.

The crude adiabatic construction is useful pedagogically, but its accuracy depends strongly on the chosen reference geometry $\(R_0\)$. Far from this geometry, the fixed electronic basis may no longer represent the relevant electronic character efficiently. More sophisticated quasi-diabatic or diabatization procedures are therefore often needed for broad regions of nuclear configuration space.

For the simulations developed in this project, the adiabatic representation is the most natural starting point because standard electronic-structure methods directly provide adiabatic potential energy surfaces, gradients, spin-orbit couplings, and nonadiabatic derivative couplings at fixed nuclear geometries. The diabatic representation remains conceptually important because it clarifies how nonadiabatic effects can alternatively be represented as potential couplings rather than derivative couplings.
