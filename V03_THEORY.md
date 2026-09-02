# v0.3 Theory: Electronic-Structure Providers and the Bridge to Direct-Dynamics Gaussian Methods

Version 0.3 does not pretend that an analytic avoided-crossing Hamiltonian is an
ab initio AIMS calculation. Instead, it introduces the missing architectural object:
a rigorously defined **electronic-structure provider**.

---

## 1. What the dynamics algorithm actually needs

At a nuclear geometry $R$, adiabatic Gaussian/nonadiabatic dynamics may require

$$
\boxed{
E_I(R),
\qquad
\nabla_R E_I(R),
\qquad
d_{IJ}(R)
=
\langle\Phi_I|\nabla_R\Phi_J\rangle.
}
$$

The propagation algorithm should depend on these physical quantities, not on the
details of the electronic-structure package that produced them.

v0.3 therefore defines a provider contract.

---

## 2. Reduction from Cartesian geometry to a pedagogical coordinate

The current framework remains one dimensional, with generalized coordinate $q$.

Let the molecular Cartesian geometry be

$$
R=R(q).
$$

A Cartesian energy gradient is

$$
\nabla_R E_I.
$$

By the chain rule,

$$
\boxed{
\frac{dE_I}{dq}
=
\nabla_R E_I
\cdot
\frac{dR}{dq}.
}
$$

Likewise, a Cartesian derivative-coupling vector obeys

$$
\boxed{
d_{IJ}^{(q)}
=
d_{IJ}^{(R)}
\cdot
\frac{dR}{dq}.
}
$$

This projection is explicitly implemented and tested.

The simplification is therefore transparent: v0.3 is not discarding the molecular
gradient/NAC information arbitrarily; it projects it onto a chosen one-dimensional
reaction coordinate.

---

## 3. Why state-averaged multireference electronic structure is relevant

Near an avoided crossing or conical-intersection region, several electronic states
must be described on a comparable footing.

State-specific orbital optimization can change character as roots approach or reorder.
State-averaged CASSCF instead optimizes a common orbital set using a weighted average

$$
\boxed{
E_{\rm SA}
=
\sum_I w_I E_I,
\qquad
\sum_Iw_I=1.
}
$$

This is why the optional PySCF provider uses an equal-weight SA-CASSCF calculation for
the requested states.

The dynamics still uses the individual state energies, state gradients, and
interstate derivative couplings—not the averaged energy as a physical surface.

---

## 4. Derivative-coupling convention

v0.3 adopts one explicit convention everywhere:

$$
\boxed{
d_{IJ}
=
\langle\Phi_I|\nabla\Phi_J\rangle.
}
$$

For real states,

$$
d_{JI}=-d_{IJ}.
$$

The optional PySCF adapter documents its conversion carefully because the PySCF
SA-CASSCF NAC interface defines `state=(ket, bra)` and returns

$$
\langle {\rm bra}|\nabla\,{\rm ket}\rangle.
$$

Therefore v0.3 requests the tuple `(J,I)` when constructing $d_{IJ}$.

This convention is recorded in provider metadata so downstream dynamics does not have
to guess.

---

## 5. Energy-gap-scaled NACs are not the same quantity

Near degeneracy,

$$
d_{IJ}
=
\frac{
\langle I|\nabla H|J\rangle
}{
E_J-E_I
}
$$

can become numerically large.

Some electronic-structure interfaces can instead return

$$
\boxed{
(E_J-E_I)d_{IJ}
=
\langle I|\nabla H|J\rangle.
}
$$

This can be numerically better behaved near a crossing, but it is **not** the
derivative coupling itself.

The v0.3 provider contract stores unscaled $d_{IJ}$. If a backend uses a scaled
quantity internally, it must convert or explicitly expose a different field rather
than silently changing semantics.

---

## 6. State tracking and gauge are provider responsibilities

Interpolation is only meaningful if the underlying data have consistent state labels
and gauge.

A tabulated provider therefore assumes that upstream data have already been:

1. state tracked;
2. phase/sign aligned;
3. placed in one declared NAC convention.

It does not silently sort states by a heuristic observable or flip signs to make plots
look smooth.

That separation is deliberate. State tracking is a scientific operation and must be
visible in the provenance.

---

## 7. Caching is part of direct-dynamics reproducibility

On-the-fly electronic-structure evaluations are expensive.

A direct-dynamics workflow should identify an electronic-structure point by at least

$$
\{
R,\ {\rm method},\ {\rm basis},\ {\rm active\ space},\ {\rm state\ averaging},
\ {\rm convergence\ settings}
\}.
$$

If the same point is requested again under the same contract, recomputing it is
unnecessary and can introduce avoidable nondeterminism.

v0.3 includes a deterministic cache wrapper and a SHA-256 point fingerprint.

A production implementation should persist the cache together with method metadata,
code versions, and input hashes.

---

## 8. Provider-guided trajectory basis functions

A TBF center on state $I$ obeys

$$
\dot q=\frac{p}{M},
$$

$$
\dot p=-\frac{dE_I}{dq}.
$$

With velocity Verlet,

$$
p_{n+1/2}
=
p_n-\frac{\Delta t}{2}E_I'(q_n),
$$

$$
q_{n+1}
=
q_n+\Delta t\frac{p_{n+1/2}}{M},
$$

$$
p_{n+1}
=
p_{n+1/2}
-\frac{\Delta t}{2}E_I'(q_{n+1}).
$$

The only difference from the analytic v0.2 trajectory guidance is where
$E_I'$ comes from: the provider rather than a hard-coded model.

---

## 9. Provider-guided spawning indicator

For a parent TBF,

$$
\boxed{
\eta_{IJ}
=
\left|
\dot q\,d_{IJ}^{(q)}
\right|
=
\left|
\frac{p}{M}d_{IJ}^{(q)}
\right|.
}
$$

The provider supplies $d_{IJ}^{(q)}$.

The same local energy-conserving child rule from v0.2 can then be used as a transparent
prototype:

$$
p_J^2
=
p_I^2+2M(E_I-E_J).
$$

Again, this is an instructional spawning rule, not the full AIMS spawning optimization.

---

## 10. What makes this "AIMS-style architecture" but not yet full AIMS

The architecture now has the correct separation:

```text
Gaussian dynamics
      |
      v
ElectronicStructureProvider.evaluate(R)
      |
      +--> energies
      +--> gradients
      +--> derivative couplings
      +--> provenance/convergence metadata
```

That is the software boundary needed for direct dynamics.

However, full AIMS additionally requires:

- multidimensional TBFs;
- robust spawned-basis placement/optimization;
- matrix-element approximations appropriate to on-the-fly electronic structure;
- coupled coefficient propagation for the dynamically changing basis;
- electronic-state phase tracking across arbitrary molecular geometries;
- initial-condition sampling;
- spawning/death/overlap controls;
- convergence with respect to all of the above.

v0.3 intentionally stops at the provider/direct-dynamics bridge rather than
mislabeling an incomplete implementation.

---

## 11. Acceptance criteria

- Analytic provider reproduces v0.2 energies, gradients, and NACs.
- Tabulated interpolation reproduces source points exactly.
- NAC antisymmetry is maintained after interpolation.
- Cartesian-to-$q$ projection obeys the chain rule on a known vector.
- Cache reuse is deterministic and measurable.
- Provider-guided velocity Verlet uses only provider-returned quantities.
- Child momentum conserves local classical energy when a real solution exists.
- PySCF is optional: absence produces a clear installation error, not a cryptic import
  failure.
- Every provider point can carry method/convergence metadata and a deterministic
  fingerprint.
