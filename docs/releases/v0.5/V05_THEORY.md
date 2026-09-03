# v0.5 Theory: Backend-Driven Multidimensional Gaussian Direct Dynamics

Version 0.5 connects the Gaussian dynamics developed in v0.1-v0.4 to an explicit
electronic-structure backend.

The central design change is:

```text
Gaussian dynamics
        |
        v
Generalized-coordinate electronic provider
        |
        v
Cartesian electronic-structure backend
        |
        +-- analytic/reference backend
        +-- PySCF SA-CASSCF backend
```

The dynamics no longer needs to know how an electronic-structure package obtains
energies, gradients, or nonadiabatic couplings.

Atomic units are used in the dynamical equations.

---

# 1. Cartesian molecular geometry

For a molecule with $N$ atoms, collect Cartesian nuclear coordinates into

$$
\boxed{
\mathbf R
=
(R_{1x},R_{1y},R_{1z},\ldots,R_{Nz})^T
\in\mathbb{R}^{3N}.
}
$$

The electronic backend evaluates, for each adiabatic state $I$,

$$
\boxed{
E_I(\mathbf R),
}
$$

the Cartesian gradient

$$
\boxed{
\mathbf g_I^{(R)}
=
\nabla_{\mathbf R}E_I,
}
$$

and derivative-coupling vectors

$$
\boxed{
\mathbf d_{IJ}^{(R)}
=
\langle\Phi_I|
\nabla_{\mathbf R}
\Phi_J\rangle.
}
$$

For real adiabatic states,

$$
\boxed{
\mathbf d_{JI}^{(R)}
=
-\mathbf d_{IJ}^{(R)}.
}
$$

The v0.5 backend contract stores these quantities explicitly with shapes

```text
energies        (nstate,)
gradients_cart  (nstate, natom, 3)
nac_cart        (nstate, nstate, natom, 3)
```

rather than hiding the coordinate convention in downstream code.

---

# 2. Generalized coordinates

A direct-dynamics Gaussian basis need not use raw Cartesian coordinates.

Let

$$
\boxed{
\mathbf R(\mathbf q)
=
\mathbf R_0 + J\mathbf q,
}
$$

where

$$
\mathbf q\in\mathbb{R}^D
$$

and

$$
J
=
\frac{\partial\mathbf R}{\partial\mathbf q}
$$

is a constant $3N\times D$ Jacobian for the linear-coordinate map implemented in
v0.5.

The chain rule gives

$$
\boxed{
\nabla_{\mathbf q}E_I
=
J^T\nabla_{\mathbf R}E_I.
}
$$

Likewise,

$$
\boxed{
\mathbf d_{IJ}^{(q)}
=
J^T\mathbf d_{IJ}^{(R)}.
}
$$

No new electronic-structure calculation is required for the projection.

---

# 3. Generalized-coordinate mass matrix

The Cartesian nuclear kinetic energy is

$$
T_N
=
\frac12
\dot{\mathbf R}^T
M_R
\dot{\mathbf R},
$$

where

$$
M_R
=
\operatorname{diag}
(m_1,m_1,m_1,\ldots,m_N,m_N,m_N).
$$

Because

$$
\dot{\mathbf R}=J\dot{\mathbf q},
$$

we obtain

$$
T_N
=
\frac12
\dot{\mathbf q}^T
J^TM_RJ
\dot{\mathbf q}.
$$

Define

$$
\boxed{
M_q=J^TM_RJ.
}
$$

The canonical momentum is

$$
\boxed{
\mathbf p=M_q\dot{\mathbf q},
}
$$

so

$$
\boxed{
\dot{\mathbf q}=M_q^{-1}\mathbf p.
}
$$

The classical kinetic energy becomes

$$
\boxed{
T_N
=
\frac12
\mathbf p^TM_q^{-1}\mathbf p.
}
$$

This formulation supports Cartesian coordinates, normal-mode-like coordinates, and
simple reaction-coordinate subspaces with one common set of equations.

---

# 4. Atomic masses and atomic units

Electronic-structure programs commonly report atomic masses in unified atomic mass
units, $u$.

The dynamical mass matrix must be expressed in electron masses when Hartree atomic
units are used.

v0.5 uses the conversion

$$
\boxed{
1\,u
\approx
1822.888486209\,m_e.
}
$$

The value is stored in one module-level constant, not duplicated throughout the code.

---

# 5. State-averaged CASSCF backend

Near an avoided crossing or conical-intersection region, independent state-specific
orbital optimizations may change character or root order.

The explicit PySCF backend therefore uses state-averaged CASSCF:

$$
\boxed{
E_{\mathrm{SA}}
=
\sum_I w_I E_I,
\qquad
w_I\ge0,
\qquad
\sum_Iw_I=1.
}
$$

The average energy is the orbital-optimization objective.

The dynamics does **not** propagate on $E_{\mathrm{SA}}$.

Instead, after convergence the backend extracts the individual

$$
E_I,
\qquad
\nabla E_I,
\qquad
\mathbf d_{IJ}.
$$

The implementation uses PySCF's `state_average_` CASSCF interface and
state-resolved SA-CASSCF gradient/NAC interfaces.

---

# 6. Explicit PySCF NAC convention

v0.5 adopts

$$
\boxed{
\mathbf d_{IJ}
=
\langle\Phi_I|
\nabla
\Phi_J\rangle.
}
$$

PySCF's SA-CASSCF NAC API defines its tuple as

```text
state = (ket, bra)
```

and returns

$$
\langle\mathrm{bra}|
\nabla
\mathrm{ket}\rangle.
$$

Therefore, to obtain the v0.5 quantity

$$
\mathbf d_{IJ}
=
\langle I|\nabla J\rangle,
$$

the backend requests

```python
state=(J, I)
```

and then explicitly writes

$$
\boxed{
\mathbf d_{JI}=-\mathbf d_{IJ}.
}
$$

This convention conversion is implemented in one location and recorded in metadata.

---

# 7. `mult_ediff` is not silently used

Near a degeneracy, the derivative coupling can become large because

$$
\mathbf d_{IJ}
=
\frac{
\langle I|\nabla H_e|J\rangle
}{
E_J-E_I
}.
$$

PySCF can return a gap-scaled NAC using `mult_ediff=True`.

That object is useful numerically, but it is not the same array required by the
dynamics equation

$$
\dot c_I
=
-iE_Ic_I
-
\sum_J
\dot{\mathbf R}\cdot
\mathbf d_{IJ}
c_J.
$$

The v0.5 dynamics therefore always requests

```python
mult_ediff=False
```

for its physical NAC field.

An optional second backend call can store PySCF's scaled quantity for diagnostics,
with its exact PySCF convention declared in metadata.

---

# 8. Electron translation factors

PySCF exposes a `use_etfs` option for SA-CASSCF NACs.

This is not treated as an invisible numerical flag.

The backend stores the selected value in every point's metadata because it changes the
definition of the returned NAC by controlling whether the relevant translational
contribution is retained.

The default in v0.5 is explicit and configurable rather than hard-coded into the
dynamics layer.

---

# 9. Electronic-structure convergence is part of the data contract

A direct-dynamics step is not valid merely because an array was returned.

The backend checks:

1. SCF convergence;
2. SA-CASSCF convergence;
3. finite energies;
4. finite gradients;
5. finite NACs;
6. NAC antisymmetry after convention conversion.

Metadata records the principal calculation settings:

- basis;
- charge and spin;
- SCF reference;
- active-space size;
- active electrons;
- number of averaged states;
- state weights;
- SCF tolerance;
- CASSCF energy tolerance;
- CASSCF orbital-gradient tolerance;
- macro-iteration limit;
- ETF choice;
- PySCF version.

A failed electronic-structure calculation raises an error rather than being converted
into a dynamics point.

---

# 10. Warm starts

A direct trajectory evaluates nearby geometries repeatedly.

PySCF permits a CASSCF calculation to start from supplied molecular orbitals.

v0.5 therefore supports an optional warm start:

$$
\boxed{
C_{\mathrm{MO}}^{(n)}
\rightarrow
\text{initial orbitals for geometry }n+1.
}
$$

This may reduce orbital-optimization cost.

However,

$$
\boxed{
\text{warm start}\neq\text{state tracking}.
}
$$

Root identity and electronic gauge remain separate scientific problems.

v0.5 does not claim that energy ordering plus warm-start orbitals gives robust global
state tracking through an arbitrary conical-intersection seam.

---

# 11. Gaussian basis in generalized coordinates

For a $D$-dimensional generalized coordinate $\mathbf q$, use a frozen Gaussian

$$
g_i(\mathbf q)
=
N
\exp
\left[
-\frac12
(\mathbf q-\mathbf q_i)^TA
(\mathbf q-\mathbf q_i)
+
i\mathbf p_i^T
(\mathbf q-\mathbf q_i)
\right].
$$

For the direct-dynamics local matrix-element layer in v0.5, interacting Gaussian pairs
are required to share the same real positive-definite width matrix $A$.

This restriction keeps the analytic matrix elements transparent.

---

# 12. Equal-width Gaussian overlap

For two Gaussians,

$$
\Delta q=q_i-q_j,
\qquad
\Delta p=p_i-p_j,
$$

the overlap is

$$
\boxed{
S_{ij}
=
\exp
\left[
-\frac14\Delta q^TA\Delta q
-\frac14\Delta p^TA^{-1}\Delta p
+
\frac{i}{2}
(p_i+p_j)^T
(q_i-q_j)
\right].
}
$$

The complex overlap centroid is

$$
\boxed{
\mu_{ij}
=
\frac{q_i+q_j}{2}
+
\frac{i}{2}
A^{-1}(p_j-p_i).
}
$$

This centroid appears naturally in all polynomial Gaussian matrix elements.

---

# 13. Gradient matrix element

Because

$$
\nabla g_j
=
[-A(q-q_j)+ip_j]g_j,
$$

the cross matrix element is

$$
\boxed{
\mathbf G_{ij}
=
\langle g_i|\nabla g_j\rangle
=
\left[
-A(\mu_{ij}-q_j)
+
ip_j
\right]
S_{ij}.
}
$$

Integration by parts gives the useful identity

$$
\boxed{
\mathbf G_{ji}^*
=
-\mathbf G_{ij}.
}
$$

That identity is what makes the first-order NAC Hamiltonian Hermitian.

---

# 14. Kinetic matrix element for a general mass matrix

Let

$$
B=M_q^{-1}.
$$

The nuclear kinetic operator is

$$
\boxed{
\hat T
=
-\frac12
\nabla^TB\nabla.
}
$$

Using integration by parts,

$$
T_{ij}
=
\frac12
\langle\nabla g_i|
B|
\nabla g_j\rangle.
$$

Define

$$
u_i
=
-A(\mu_{ij}-q_i)-ip_i,
$$

$$
u_j
=
-A(\mu_{ij}-q_j)+ip_j.
$$

The cross-Gaussian covariance is

$$
\frac12A^{-1}.
$$

Therefore,

$$
\boxed{
T_{ij}
=
\frac12
S_{ij}
\left[
u_i^TBu_j
+
\frac12
\operatorname{Tr}(BA)
\right].
}
$$

For $i=j$,

$$
\boxed{
T_{ii}
=
\frac12p_i^TBp_i
+
\frac14\operatorname{Tr}(BA),
}
$$

which is the classical center kinetic energy plus the Gaussian zero-point kinetic
contribution.

---

# 15. Local electronic approximation

A true multidimensional AIMS matrix element requires electronic quantities throughout
the nuclear overlap region.

That is too expensive to evaluate exactly on the fly.

v0.5 therefore introduces a deliberately named approximation:

> **local constant-electronic-quantity Gaussian matrix approximation**

For a Gaussian pair $i,j$, evaluate the electronic provider once at

$$
\boxed{
q_{ij}
=
\frac{q_i+q_j}{2}.
}
$$

Over the overlap region, approximate

$$
E_a(q)\approx E_a(q_{ij}),
$$

and

$$
d_{ab}(q)\approx d_{ab}(q_{ij}).
$$

The approximation is visible and independently testable.

---

# 16. Locally constant covariant kinetic operator

The adiabatic kinetic operator can be written

$$
\boxed{
\hat T_{\mathrm{ad}}
=
-\frac12
(\nabla+d)^T
B
(\nabla+d).
}
$$

If $d$ is treated as locally constant, its divergence is neglected over the Gaussian
pair overlap region.

Expanding gives

$$
\hat T_{\mathrm{ad}}
\approx
-\frac12\nabla^TB\nabla
-
d^TB\nabla
-
\frac12d^TBd.
$$

Here $d$ is a matrix in electronic-state space and a vector in nuclear-coordinate
space.

Define

$$
\boxed{
D^{(2)}_{ab}
=
\sum_{c,\alpha,\beta}
d_{ac,\alpha}
B_{\alpha\beta}
d_{cb,\beta}.
}
$$

Then the local Gaussian Hamiltonian matrix is

$$
\boxed{
H_{ia,jb}
=
\delta_{ab}
\left[
T_{ij}
+
E_a(q_{ij})S_{ij}
\right]
-
\sum_{\alpha\beta}
d_{ab,\alpha}
B_{\alpha\beta}
G_{ij,\beta}
-
\frac12
D^{(2)}_{ab}
S_{ij}.
}
$$

This includes:

- exact equal-width Gaussian kinetic integrals;
- centroid adiabatic energies;
- first-order derivative-coupling terms;
- the locally constant $d^2$ contribution.

It neglects

$$
\nabla\cdot(Bd),
$$

because the electronic coupling is frozen at the pair centroid.

That omission is part of the declared approximation, not an accidental missing term.

---

# 17. Hermiticity of the local Hamiltonian

At a common pair centroid,

$$
d_{ba}=-d_{ab},
$$

and

$$
G_{ji}^*=-G_{ij}.
$$

Therefore the first-order NAC terms satisfy

$$
H_{jb,ia}^*
=
H_{ia,jb}.
$$

The $d^2$ matrix is Hermitian for a real antisymmetric derivative-coupling field and
real symmetric mass metric.

Consequently,

$$
\boxed{
H=H^\dagger
}
$$

within the local approximation.

The test suite verifies this numerically for synthetic provider data.

---

# 18. Moving-basis matrix

For a fixed-width Gaussian,

$$
\dot g_j
=
\left[
(A(q-q_j)-ip_j)^T\dot q_j
+
i(q-q_j)^T\dot p_j
\right]g_j.
$$

Using the complex centroid,

$$
\boxed{
T^{\mathrm{basis}}_{ij}
=
\langle g_i|\dot g_j\rangle
}
$$

becomes

$$
\boxed{
T^{\mathrm{basis}}_{ij}
=
S_{ij}
\left[
(A(\mu_{ij}-q_j)-ip_j)^T\dot q_j
+
i(\mu_{ij}-q_j)^T\dot p_j
\right].
}
$$

The electronic overlap makes this term zero between different adiabatic state labels
in the local TBF approximation.

As before,

$$
\boxed{
\dot S
=
T^{\mathrm{basis}}
+
T^{\mathrm{basis}\dagger}.
}
$$

---

# 19. Coupled coefficient equation

The dynamically changing nonorthogonal Gaussian basis obeys

$$
\boxed{
iS\dot C
=
(H-iT^{\mathrm{basis}})C.
}
$$

The code solves the linear system

$$
S\dot C
=
-iHC
-
T^{\mathrm{basis}}C
$$

rather than constructing $S^{-1}$ explicitly.

The physical norm is

$$
\boxed{
\|\Psi\|^2
=
C^\dagger SC.
}
$$

---

# 20. Backend-driven classical guidance

For a TBF on state $a$,

$$
\boxed{
\dot q
=
M_q^{-1}p,
}
$$

$$
\boxed{
\dot p
=
-\nabla_q E_a.
}
$$

The gradient comes directly from the projected electronic-structure backend.

Thus the same TBF propagation code can use:

- an analytic benchmark provider;
- tabulated data;
- a cached PySCF provider.

---

# 21. Backend-driven spawning indicator

The nonadiabatic scalar along the TBF velocity is

$$
\boxed{
\eta_{ab}
=
|
\dot q^T d_{ab}
|.
}
$$

Since

$$
\dot q=M_q^{-1}p,
$$

$$
\boxed{
\eta_{ab}
=
|
p^TM_q^{-1}d_{ab}
|.
}
$$

This is the multidimensional generalization of the 1D and 2D indicators used in
earlier versions.

---

# 22. Energy-conserving child momentum with a mass metric

Let the child momentum be

$$
p_b=p_a+\lambda n,
$$

where $n$ is the Euclidean-normalized derivative-coupling direction.

Require

$$
\frac12
p_b^TBp_b
+
E_b
=
\frac12
p_a^TBp_a
+
E_a.
$$

Define

$$
A_n=n^TBn,
$$

$$
B_n=p_a^TBn,
$$

and

$$
\Delta E=E_b-E_a.
$$

Then

$$
\boxed{
A_n\lambda^2
+
2B_n\lambda
+
2\Delta E
=
0.
}
$$

The discriminant is

$$
\boxed{
\mathcal D
=
B_n^2
-
2A_n\Delta E.
}
$$

If

$$
\mathcal D<0,
$$

no real locally energy-conserving child exists under this placement rule.

Otherwise,

$$
\boxed{
\lambda
=
\frac{
-B_n\pm\sqrt{\mathcal D}
}{
A_n
}.
}
$$

v0.5 chooses the real root producing the smaller $|\lambda|$.

---

# 23. Direct-dynamics cache

An ab initio Gaussian calculation repeatedly requests electronic data at TBF centers
and Gaussian-pair centroids.

A cache key therefore contains:

- provider namespace/settings fingerprint;
- generalized coordinate;
- deterministic coordinate rounding.

A cached result stores

$$
E,
\quad
\nabla_qE,
\quad
d_q,
\quad
M_q,
\quad
\text{metadata}.
$$

The cache does not change the physics.

It makes repeated matrix construction reproducible and prevents identical
electronic-structure points from being recomputed unnecessarily.

---

# 24. What PySCF v0.5 does explicitly

The explicit backend executes the following stages:

```text
MolecularGeometry in bohr
        |
        v
gto.M(..., unit="Bohr")
        |
        v
RHF or ROHF
        |
        v
CASSCF(active space)
        |
        v
state_average_(weights)
        |
        v
SA-CASSCF kernel
        |
        +--> individual state energies
        +--> state-specific SA-CASSCF gradients
        +--> pairwise SA-CASSCF NACs
        +--> atomic masses
        +--> convergence + method metadata
```

No placeholder electronic energies are inserted after PySCF is selected.

---

# 25. What v0.5 still does not solve

The following remain beyond the current release:

1. rigorous cross-geometry electronic-state tracking from many-electron wavefunction
   overlaps;
2. multidimensional nonlinear internal-coordinate metric tensors;
3. changing Gaussian width matrices during the direct spawned-basis run;
4. production AIMS optimal spawning;
5. full saddle-point/SPA0/SPA1 hierarchy;
6. adaptive electronic/nuclear time stepping;
7. large initial-condition ensembles;
8. production-scale pruning and regularization;
9. robust dynamics exactly on a CI seam where individual adiabatic states are not
   uniquely defined.

The purpose of v0.5 is narrower and concrete:

$$
\boxed{
\text{replace the analytic electronic model by a real, explicit,
validated electronic-structure backend boundary.}
}
$$
