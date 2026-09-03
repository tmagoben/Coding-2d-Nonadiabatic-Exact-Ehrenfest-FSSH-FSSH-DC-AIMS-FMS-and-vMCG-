# v0.14 Theory: Fully Time-Adaptive TDSE-Defect Control

Version 0.15 deliberately folds the planned v0.14 control-loop work into the next
release rather than creating an artificial intermediate package.

The progression is therefore

```text
v0.13
static initial residual selection
+ instantaneous defect-enrichment primitive

        ↓

v0.14
fully time-adaptive defect checks
+ residual-triggered growth
+ low-loss pruning
+ hysteresis/cooldowns
+ explicit algorithmic-complexity accounting
+ Hermitian half-build optimization
```

The central goal is to turn the v0.13 residual mathematics into an actual propagation
controller.

Atomic units are used throughout.

---

## 1. Adaptive state

At time $t$, the spinor-complete Gaussian approximation is

$$
\boxed{
\Psi_G(R,t)
=
\sum_{i=1}^{N(t)}
g_i(R,t)\mathbf C_i(t).
}
$$

The basis dimension

$$
N(t)
$$

is now itself time dependent.

Each Gaussian carries a complete two-state electronic coefficient vector, while its
integer `state` is used only for classical guidance of the Gaussian center.

---

## 2. Propagation between adaptation events

Between basis changes,

$$
iS\dot C=(H-iT)C.
$$

The midpoint/Cayley update remains

$$
\boxed{
\left[
S_m+\frac{\Delta t}{2}(iH_m+T_m)
\right]C_{n+1}
=
\left[
S_m-\frac{\Delta t}{2}(iH_m+T_m)
\right]C_n.
}
$$

Thus v0.14 does not replace the validated v0.12-v0.13 coefficient propagation.

It adds a control layer around it.

---

## 3. Error monitor

Every

$$
m
$$

steps, v0.14 reconstructs

$$
\boxed{
\mathcal R
=
i\dot\Psi_G-H\Psi_G
}
$$

on the independent diagnostic grid.

The principal normalized trigger is

$$
\boxed{
\eta
=
\frac{\|\mathcal R\|}{\|H\Psi_G\|}.
}
$$

This is the quantity called `relative_to_hpsi` in the code.

---

## 4. Hysteresis

A single threshold would cause unstable basis-size chatter:

```text
add a TBF
defect falls slightly
remove a TBF
defect rises
add it again
...
```

v0.14 therefore uses two thresholds:

$$
\eta_{\mathrm{add}}
>
\eta_{\mathrm{remove}}.
$$

For the release benchmark,

$$
\boxed{
\eta_{\mathrm{add}}=0.020,
\qquad
\eta_{\mathrm{remove}}=0.006.
}
$$

The interval

$$
0.006<\eta<0.020
$$

is a no-action band.

This is classical control-system hysteresis applied to basis adaptation.

---

## 5. Adaptation cooldown

After a growth or pruning event, the controller waits a minimum number of propagation
steps before another ordinary adaptation.

For the release,

$$
\boxed{
\Delta n_{\mathrm{adapt,min}}=10\;\text{steps}.
}
$$

This separates the physical response to one basis change from the decision to make a
second change.

---

## 6. Growth criterion

If

$$
\eta\ge\eta_{\mathrm{add}},
$$

the code constructs a local candidate dictionary around the current TBFs.

Candidates are generated with the existing energy-conserving local placement
machinery, but the old coupling-based rank is discarded.

For candidate $c$,

$$
g_c^\perp=(I-P_\mathcal B)g_c
$$

and the TDSE-defect capture is

$$
\boxed{
\Delta_c^{\mathrm{TDSE}}
=
\frac{
\sum_a
|\langle g_c^\perp|\mathcal R_a\rangle|^2
}{\|g_c^\perp\|^2}.
}
$$

The candidate with the largest admissible capture fraction is selected.

---

## 7. Physically guided candidate dictionary

v0.14 separates **candidate generation** from **candidate ranking**.

Generation uses physically interpretable constraints:

- parent TBF center;
- same- or other-surface future guidance;
- small local position shifts;
- width changes;
- momentum adjustment along NAC or momentum directions;
- local classical-energy conservation.

Ranking uses the measured TDSE defect.

Therefore the algorithm is neither

```text
pure coupling-triggered spawning
```

nor

```text
unconstrained global Gaussian optimization.
```

It is a physically constrained residual search.

---

## 8. Zero-amplitude insertion

An accepted Gaussian enters with

$$
\boxed{
C_{c0}=C_{c1}=0.
}
$$

Therefore

$$
\Psi_{\mathrm{after}}(t_*)
=
\Psi_{\mathrm{before}}(t_*).
$$

The basis change alters only the future Galerkin tangent space.

The release event at step 10 reduced the normalized defect from

$$
0.032380191
$$

to

$$
\boxed{
0.029987958
}.
$$

The predicted defect capture fraction was

$$
\boxed{
0.1506769.
}
$$

---

## 9. Exact low-loss pruning

v0.14 adds a new pruning score that is exact for deleting one nuclear Gaussian from
the represented spinor-complete wavefunction.

Let

$$
S
$$

be the nuclear overlap matrix and

$$
\mathbf C_j
$$

the complete electronic coefficient row on Gaussian $j$.

The component of $g_j$ orthogonal to all other Gaussians has norm

$$
\boxed{
n_j=
\frac{1}{(S^{-1})_{jj}}.
}
$$

The exact squared represented-wavefunction loss from optimally deleting $j$ is

$$
\boxed{
L_j
=
n_j
\sum_a|C_{ja}|^2.
}
$$

All $N$ deletion scores are obtained from one dense solve/inverse of $S$.

---

## 10. Pruning policy

Pruning can occur for three different reasons:

1. **sustained low defect** — the basis is larger than needed;
2. **basis-budget replacement** — the defect is high but the hard basis cap has been
   reached;
3. **emergency conditioning** — the overlap matrix exceeds the hard condition limit.

Recently added Gaussians can be protected for a minimum age to prevent immediate
add/remove oscillation.

---

## 11. Pruning stress validation

The release deliberately adds a nearly redundant, zero-amplitude Gaussian and then
audits the deletion rule.

The result is:

```text
fractional wavefunction loss:
0.0

condition number before:
379319.12346481933

condition number after:
67.30166373596352

condition improvement factor:
5636.103216600353
```

Thus the stress-test removal changes the represented wavefunction by zero while
improving conditioning by more than three orders of magnitude.

---

## 12. Hermitian half-build

For the exact LVC Gaussian matrices,

$$
S_{ji}=S_{ij}^*
$$

and

$$
H_{ji}=H_{ij}^\dagger.
$$

Earlier builders evaluated all ordered Gaussian pairs:

$$
N^2.
$$

v0.14 evaluates only

$$
\boxed{
\frac{N(N+1)}{2}.
}
$$

The fractional pair-evaluation saving is

$$
\boxed{
1-
\frac{N(N+1)/2}{N^2}
=
\frac12-\frac{1}{2N}.
}
$$

At

$$
N=11,
$$

this is

$$
45.45\%.
$$

The actual adaptive campaign reduced pair evaluations by

$$
\boxed{
100\times0.4542\%.
}
$$

The asymptotic class remains

$$
O(N^2),
$$

but the leading pair-work constant is almost halved.

---

## 13. Release benchmark

The v0.14 benchmark starts from a 10-Gaussian residual-selected initial bank.

It grows to 11 Gaussians at

$$
t=0.05.
$$

The time-averaged basis size is

$$
\boxed{
\bar N=10.925.
}
$$

The final representation-consistent projected-state error is

$$
\boxed{
\epsilon_{\mathrm{dyn}}
=
9.5278046e-05.
}
$$

The final original-target density error is

$$
\boxed{
\epsilon_\rho
=
0.03330494.
}
$$

---

## 14. Relation to v0.13

The v0.13 reference used a static 11-Gaussian bank.

v0.14 begins with 10 and decides during propagation whether the eleventh direction is
actually needed.

The comparison is:

```text
v0.13 target density error:
0.03178630139393256

v0.14 target density error:
0.03330494031479218

v0.13 projected dynamics error:
0.00011354880287339317

v0.14 projected dynamics error:
9.527804623132635e-05
```

v0.14 is not presented as a dramatic accuracy jump.

Its advance is that the basis dimension is now controlled by a measured dynamical
error rather than fixed before propagation.

---

## 15. Why the adaptive benchmark can be slightly less accurate than a larger static bank

A static 11-Gaussian initial bank has access to all 11 directions from

$$
t=0.
$$

The v0.14 run intentionally begins with only 10.

The eleventh Gaussian is introduced at

$$
t=0.05.
$$

Therefore some representational flexibility is intentionally deferred.

This is the basic accuracy-versus-cost tradeoff of adaptive basis methods.

The correct validation question is not

> Does adaptation always beat a larger static basis in every scalar metric?

but

> Does the error controller add/remove directions predictably while respecting the
> chosen accuracy, conditioning, and computational budgets?

For the configured release tolerances, it does.

---

## 16. Computational complexity is now part of the numerical contract

Every v0.14 adaptive result stores:

- pair-matrix build count;
- actual Hermitian pair evaluations;
- ordered-pair equivalent;
- moving-basis $T$-matrix count;
- dense Cayley solve count;
- TDSE-defect evaluation count;
- candidate ranking count;
- total candidates scored;
- enrichment/pruning event counts;
- peak basis/electronic dimensions;
- timing by algorithmic category;
- symbolic Big-O scaling.

The release therefore treats computational complexity as an observable of the
algorithm rather than an undocumented implementation detail.

See `V14_ALGORITHM_COMPLEXITY.md`.

---

## 17. Scientific label

The appropriate description is:

> **time-adaptive, TDSE-defect-controlled, spinor-complete Gaussian dynamics with
> residual-driven growth, exact low-loss pruning, and explicit complexity accounting
> on an analytic LVC conical-intersection benchmark.**

It should not be called:

- production AIMS;
- full vMCG;
- a scalable full-dimensional TDSE solver;
- a completed PySCF residual-adaptation engine.
