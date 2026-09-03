# References and their role in this repository

## Heller 1975 — thawed Gaussian dynamics

E. J. Heller, **Time-dependent approach to semiclassical dynamics**,
*Journal of Chemical Physics* **62**, 1544-1555 (1975).

DOI: https://doi.org/10.1063/1.430620

Used for:

- the time-dependent Gaussian ansatz;
- local quadratic expansion of the potential;
- classical evolution of Gaussian centers;
- width and phase evolution;
- the idea that Gaussian packets provide a trajectory-like representation of
  semiclassical dynamics.

## Heller 1975 — Gaussian phase-space basis

E. J. Heller, **Wavepacket path integral formulation of semiclassical dynamics**,
*Chemical Physics Letters* **34**, 321-325 (1975).

DOI: https://doi.org/10.1016/0009-2614(75)85284-5

Used for:

- motivation for representing a more general wavefunction as a superposition of
  phase-space-localized Gaussian packets.

## Martínez, Ben-Nun, and Levine 1996 — multiple electronic states

T. J. Martínez, M. Ben-Nun, and R. D. Levine,
**Multi-Electronic-State Molecular Dynamics: A Wave Function Approach with
Applications**,
*Journal of Physical Chemistry* **100**, 7884-7895 (1996).

DOI: https://doi.org/10.1021/jp953105a

Used for:

- the wavefunction-based viewpoint underlying multiple-spawning dynamics;
- motivation for permitting distinct nuclear dynamics on different electronic
  states while retaining coherent amplitudes.

## Ben-Nun, Quenneville, and Martínez 2000 — AIMS

M. Ben-Nun, J. Quenneville, and T. J. Martínez,
**Ab Initio Multiple Spawning: Photochemistry from First Principles Quantum Molecular
Dynamics**,
*Journal of Physical Chemistry A* **104**, 5161-5175 (2000).

DOI: https://doi.org/10.1021/jp994174i

Used for:

- the AIMS framework;
- simultaneous nuclear dynamics and electronic-structure evaluation;
- adaptive frozen-Gaussian basis growth in nonadiabatic regions.

## Worth and Burghardt 2003 — variational Gaussian dynamics

G. A. Worth and I. Burghardt,
**Full quantum mechanical molecular dynamics using Gaussian wavepackets**,
*Chemical Physics Letters* **368**, 502-508 (2003).

DOI: https://doi.org/10.1016/S0009-2614(02)01920-6

Used for:

- the idea of variationally coupled Gaussian basis functions;
- the connection between moving Gaussians and converged quantum dynamics;
- motivation for the vMCG layer.

## Lasorne et al. 2006/2007 — direct vMCG

B. Lasorne, M. J. Bearpark, M. A. Robb, and G. A. Worth,
**Direct quantum dynamics using variational multi-configuration Gaussian wavepackets**,
*Chemical Physics Letters* **432**, 604-609 (2006).

DOI: https://doi.org/10.1016/j.cplett.2006.10.099

B. Lasorne, M. A. Robb, and G. A. Worth,
**Direct quantum dynamics using variational multi-configuration Gaussian wavepackets.
Implementation details and test case**,
*Physical Chemistry Chemical Physics* (2007).

Used for:

- direct/on-the-fly vMCG perspective;
- distinction between variational Gaussian motion and classically guided Gaussian
  centers.

## Richings et al. 2015 — vMCG review

G. W. Richings, I. Polyak, K. E. Spinlove, G. A. Worth, I. Burghardt, and B. Lasorne,
**Quantum dynamics simulations using Gaussian wavepackets: the vMCG method**,
*International Reviews in Physical Chemistry* **34**, 269-308 (2015).

DOI: https://doi.org/10.1080/0144235X.2015.1051354

Used as a review-level source for the mathematical and computational position of vMCG
within Gaussian quantum dynamics.

---

## Important terminology note

The repository uses:

- **moving-basis MCG foundation** for classically/externally guided frozen Gaussian
  centers with quantum coefficients;
- **variational multi-Gaussian TDVP** for the compact real-parameter McLachlan
  implementation in this repository;
- **vMCG** when discussing the established variational multi-configuration Gaussian
  methodology in the literature;
- **spawning foundation** for the minimal basis-growth demonstration;
- **FMS/AIMS** only when discussing the established full methodologies.

This distinction prevents a small educational implementation from being mislabeled as
a production implementation of a much larger research method.


# v0.4 references: conical intersections, geometric phase, and multiple spawning

## Mead and Truhlar 1979 — conical intersections and nuclear geometric phase

C. A. Mead and D. G. Truhlar,
**On the determination of Born-Oppenheimer nuclear motion wave functions including
complications due to conical intersections and identical nuclei**,
*Journal of Chemical Physics* **70**, 2284-2296 (1979).

DOI: https://doi.org/10.1063/1.437734

Used for:

- the phase complication caused by conical intersections;
- the vector-potential/gauge viewpoint for adiabatic nuclear motion;
- the requirement that electronic and nuclear phases be handled consistently.

## Ben-Nun, Quenneville, and Martínez 2000 — AIMS

M. Ben-Nun, J. Quenneville, and T. J. Martínez,
**Ab Initio Multiple Spawning: Photochemistry from First Principles Quantum Molecular
Dynamics**,
*Journal of Physical Chemistry A* **104**, 5161-5175 (2000).

DOI: https://doi.org/10.1021/jp994174i

Used for:

- the multiple-spawning wavefunction viewpoint;
- dynamically expanded Gaussian trajectory basis functions;
- the distinction between coupled quantum amplitudes and classical trajectory labels;
- the direct-dynamics connection between nuclear propagation and electronic
  structure.

## Worth and Burghardt 2003 — variational Gaussian dynamics

G. A. Worth and I. Burghardt,
**Full quantum mechanical molecular dynamics using Gaussian wavepackets**,
*Chemical Physics Letters* **368**, 502-508 (2003).

DOI: https://doi.org/10.1016/S0009-2614(02)01920-6

Used for:

- the role of coupled moving Gaussian basis functions in converged quantum dynamics;
- the contrast between variational Gaussian motion and independent classical
  trajectory guidance.

## Lasorne et al. 2006 — direct vMCG

B. Lasorne, M. J. Bearpark, M. A. Robb, and G. A. Worth,
**Direct quantum dynamics using variational multi-configuration Gaussian
wavepackets**,
*Chemical Physics Letters* **432**, 604-609 (2006).

DOI: https://doi.org/10.1016/j.cplett.2006.10.099

Used for:

- the direct-dynamics connection between moving Gaussian quantum dynamics and
  on-the-fly electronic structure.

## Modern Gaussian nonadiabatic platform context

**Legion: A Platform for Gaussian Wavepacket Nonadiabatic Dynamics**,
*Journal of Chemical Theory and Computation*.

DOI: https://doi.org/10.1021/acs.jctc.4c01697

Used only as modern context showing that Gaussian wavepacket nonadiabatic methods,
including multiple-spawning descendants and related Gaussian trajectory methods,
remain an active software-development area.


# v0.5 implementation references

## PySCF MCSCF documentation

PySCF, **Multi-configuration self-consistent field (MCSCF)**.

Official documentation:
https://pyscf.org/user/mcscf.html

Used for:

- CASSCF construction;
- `state_average_` behavior;
- state-average weights;
- CASSCF convergence controls;
- orbital initial guesses.

## PySCF analytical NAC documentation

PySCF, **pyscf.nac package — Analytical Nonadiabatic Coupling Vectors**.

Official API:
https://pyscf.org/pyscf_api_docs/pyscf.nac.html

Used for:

- SA-CASSCF NAC backend calls;
- the explicit `state=(ket,bra)` convention;
- the meaning of `mult_ediff`;
- the `use_etfs` option.

## PySCF molecular geometry API

PySCF, **Molecular structure / gto**.

Official documentation:
https://pyscf.org/user/gto.html

Used for:

- explicit Cartesian molecule construction;
- Bohr coordinate handling.

PySCF `Mole.atom_mass_list()` documentation:
https://pyscf.org/pyscf_api_docs/pyscf.gto.html

Used for atomic masses returned with each backend point.

## CODATA/NIST constants

NIST, **2022 CODATA Recommended Values of the Fundamental Physical Constants**.

https://physics.nist.gov/constants

Used as the reference family for the unified-atomic-mass/electron-mass conversion
employed by the atomic-unit nuclear mass layer.


# v0.6 implementation references: electronic-state overlaps and tracking

## PySCF FCI overlap

PySCF, `pyscf.fci.addons.overlap`.

Official API:
https://pyscf.org/pyscf_api_docs/pyscf.fci.html

Used for:

- overlap between CI wavefunctions expressed in nonorthogonal one-particle orbital
  bases;
- the core+active many-electron SA-CASSCF overlap matrix used for root tracking.

The API accepts an overlap matrix `s` for the nonorthogonal one-particle basis.

## PySCF cross-molecule AO overlaps

PySCF, `pyscf.gto.mole.intor_cross`.

Official API:
https://pyscf.org/pyscf_api_docs/pyscf.gto.html

Used for:

$$
S_{\mu\nu}^{AB}
=
\langle\chi_\mu(\mathbf R_A)|\chi_\nu(\mathbf R_B)\rangle.
$$

This is transformed into the previous/current core+active MO overlap matrix before the
many-electron CI overlap is evaluated.

## PySCF CASSCF state-averaged data structures

PySCF MCSCF documentation:
https://pyscf.org/user/mcscf.html

PySCF MCSCF API:
https://pyscf.org/pyscf_api_docs/pyscf.mcscf.html

Used for:

- state-averaged CASSCF;
- `mc.e_states`;
- `mc.ci`;
- `mc.mo_coeff`;
- `mc.ncore`;
- `mc.ncas`;
- `mc.nelecas`.

PySCF documents `nelecas` as the active `(nalpha, nbeta)` electron tuple and supports
a list of CI vectors for state-averaged calculations.

## PySCF determinant/string utilities

PySCF, `pyscf.fci.cistring`.

Official API:
https://pyscf.org/pyscf_api_docs/pyscf.fci.html

Used for:

- determinant bit strings;
- determinant addresses;
- exact embedding of active-space CI coefficients into the restricted core+active FCI
  determinant space.

## State-tracking methodology

The implementation uses a maximum-overlap principle at the **many-electron
wavefunction** level and a unitary Procrustes/polar-decomposition construction for
local subspace alignment.

These are standard linear-algebra constructions rather than a claim that one unique
state label exists at an exact degeneracy. At exact degeneracy, v0.6 treats the
electronic subspace as the physically meaningful object.

# v0.7 references: geometric/non-Abelian gauge structure and local diabatization

## Berry 1984 — geometric phase

M. V. Berry,
**Quantal phase factors accompanying adiabatic changes**,
*Proceedings of the Royal Society of London A* **392**, 45–57 (1984).

DOI: 10.1098/rspa.1984.0023

Used for:

- closed-loop geometric phase;
- gauge-invariant loop holonomy;
- the sign change of a real eigenstate around a degeneracy.

## Wilczek and Zee 1984 — non-Abelian adiabatic gauge structure

F. Wilczek and A. Zee,
**Appearance of Gauge Structure in Simple Dynamical Systems**,
*Physical Review Letters* **52**, 2111–2114 (1984).

DOI: 10.1103/PhysRevLett.52.2111

Used for:

- matrix-valued adiabatic gauge connections in degenerate manifolds;
- non-Abelian holonomy as the natural extension of the Berry phase.

## Granucci, Persico, and Toniolo 2001 — local diabatization in dynamics

G. Granucci, M. Persico, and A. Toniolo,
**Direct semiclassical simulation of photochemical processes with semiempirical wave
functions**,
*Journal of Chemical Physics* **114**, 10608–10615 (2001).

DOI: 10.1063/1.1376633

Used for contextual motivation that wavefunction-overlap/local-diabatic electronic
propagation is a practical alternative to directly integrating sharply peaked
adiabatic derivative couplings along a trajectory.

## Ben-Nun, Quenneville, and Martínez 2000 — AIMS

M. Ben-Nun, J. Quenneville, and T. J. Martínez,
**Ab Initio Multiple Spawning: Photochemistry from First Principles Quantum Molecular
Dynamics**,
*Journal of Physical Chemistry A* **104**, 5161–5175 (2000).

DOI: 10.1021/jp994174i

Used for:

- the coupled Gaussian trajectory-basis viewpoint;
- dynamic basis growth;
- the motivation for handling electronic information consistently across multiple
  simultaneously active TBFs.

# v0.8 implementation context

v0.8 continues to use the primary methodological sources already cited for AIMS/FMS, variational Gaussian dynamics, conical-intersection geometric phase, PySCF SA-CASSCF/NACs, and nonorthogonal many-electron overlaps.

The new overlap/local-diabatic propagator is derived directly in `docs/releases/v0.8/V08_THEORY.md` from the finite-step electronic overlap relation

$$
O(t,t+\Delta t)=I+\Delta t\,\dot R\cdot d+\mathcal O(\Delta t^2),
$$

and is validated numerically against the explicit derivative-coupling equation on the repository's analytic two-state CI model.

# v0.9 references: AIMS approximations and convergence

## Curchod, Glover, and Martínez 2020 — SSAIMS / practical AIMS approximations

B. F. E. Curchod, W. J. Glover, and T. J. Martínez,
**SSAIMS—Stochastic-Selection Ab Initio Multiple Spawning for Efficient Nonadiabatic
Molecular Dynamics**, *J. Phys. Chem. A* **124**, 6133–6143 (2020).

DOI: https://doi.org/10.1021/acs.jpca.0c04113

Used for the documented AIMS context in which coupled TBFs are propagated and the
required integrals are commonly simplified with a zeroth-order saddle-point
approximation evaluated at TBF-pair centroids.

## Ibele and Curchod 2021 — SPA order and CI Hamiltonian consistency

L. M. Ibele and B. F. E. Curchod,
**Dynamics near a conical intersection—A diabolical compromise for the approximations
of ab initio multiple spawning**, *J. Chem. Phys.* **155**, 174119 (2021).

DOI: https://doi.org/10.1063/5.0071376

Used for the caution that improving AIMS matrix elements near conical intersections is
not equivalent to independently adding isolated higher-order terms: SPA order,
nonadiabatic terms, geometric-phase-related contributions, and Hamiltonian Hermiticity
must be considered consistently.

## Legion Gaussian-wavepacket platform

**Legion: A Platform for Gaussian Wavepacket Nonadiabatic Dynamics**,
*J. Chem. Theory Comput.*

DOI: https://doi.org/10.1021/acs.jctc.4c01697

Used as modern software context for AIMS/FMS, spawning/basis management, and efficient
Gaussian nonadiabatic dynamics implementations.


# v0.11 references: basis completeness and spawning

## Ben-Nun, Quenneville, and Martínez 2000 — AIMS

M. Ben-Nun, J. Quenneville, and T. J. Martínez,
**Ab Initio Multiple Spawning: Photochemistry from First Principles Quantum Molecular
Dynamics**,
*Journal of Physical Chemistry A* **104**, 5161-5175 (2000).

DOI: https://doi.org/10.1021/jp994174i

Used for:

- the coupled Gaussian TBF representation;
- adaptive basis expansion by spawning;
- the direct-dynamics/AIMS framework.

## Yang, Coe, Kaduk, and Martínez 2009 — optimal spawning

S. Yang, J. D. Coe, B. Kaduk, and T. J. Martínez,
**An "optimal" spawning algorithm for adaptive basis set expansion in nonadiabatic
dynamics**,
*Journal of Chemical Physics* **130**, 134113 (2009).

DOI: https://doi.org/10.1063/1.3103930

Used for the central v0.11 design principle:

- newly spawned basis functions should be chosen to maximize useful parent-child
  coupling;
- parent and child classical energies are constrained to agree;
- improved spawning can reduce the number of basis functions needed for convergence.

v0.11 does **not** claim to reproduce the complete continuous optimal-spawning
algorithm.  It implements a finite local, energy-constrained search inspired by that
principle.

## Vindel-Zandbergen et al. 2022 — multiple-spawning basis control

P. Vindel-Zandbergen et al.,
**Extending the Applicability of the Multiple-Spawning Framework for Nonadiabatic
Molecular Dynamics**,
*Journal of Physical Chemistry Letters* **13** (2022).

DOI: https://doi.org/10.1021/acs.jpclett.2c03295

Used for modern context on:

- adaptive TBF-basis growth;
- the computational cost of uncontrolled spawning;
- the need for principled TBF selection/removal.

## Mignolet et al. 2020 — stochastic-selection AIMS

B. Mignolet et al.,
**SSAIMS—Stochastic-Selection Ab Initio Multiple Spawning for Efficient Nonadiabatic
Molecular Dynamics**,
*Journal of Physical Chemistry A* **124**, 6133-6143 (2020).

DOI: https://doi.org/10.1021/acs.jpca.0c04113

Used for:

- the importance of TBF phase-space overlap and Hamiltonian coupling in determining
  whether Gaussian branches remain dynamically connected;
- basis-control context.

## Burghardt, Meyer, and co-workers — variational Gaussian comparison

G. A. Worth / I. Burghardt and related Gaussian-wavepacket work, including

**Multimode quantum dynamics using Gaussian wavepackets: The Gaussian-based
multiconfiguration time-dependent Hartree (G-MCTDH) method applied to the absorption
spectrum of pyrazine**,
*Journal of Chemical Physics* **129**, 174104 (2008).

DOI: https://doi.org/10.1063/1.2996349

Used only for contrast:

- variational Gaussian methods can evolve Gaussian parameters and can converge toward
  exact wavepacket dynamics;
- v0.11's width bank is **not** a vMCG/G-MCTDH variational width equation.

## Legion 2025 — modern Gaussian nonadiabatic software context

**Legion: A Platform for Gaussian Wavepacket Nonadiabatic Dynamics**,
*Journal of Chemical Theory and Computation*.

DOI: https://doi.org/10.1021/acs.jctc.4c01697

Used for modern implementation context on:

- AIMS spawning at strong nonadiabatic coupling;
- child energy conservation;
- momentum adjustment;
- modular separation of Gaussian propagation and basis management.


# v0.12 references: representation consistency, local diabatization, and coherence

## Gu 2023 — discrete-variable local diabatic representation

B. Gu,
**A Discrete-Variable Local Diabatic Representation of Conical Intersection
Dynamics**,
*Journal of Chemical Theory and Computation* **19**, 6557-6563 (2023).

DOI: https://doi.org/10.1021/acs.jctc.3c00560

Used as primary literature context for:

- avoiding the singular first- and second-derivative couplings of a globally
  adiabatic CI representation;
- using electronic overlap information in a local-diabatic representation;
- retaining nonadiabatic transitions, electronic coherence, and geometric phase.

The v0.12 analytic global-diabatic Gaussian benchmark is not an implementation of
Gu's DVR method.  The reference supports the broader representation principle.

## Local diabatic representation with Strang splitting and Fourier basis

B. Gu,
**Nonadiabatic Conical Intersection Dynamics in the Local Diabatic Representation
with Strang Splitting and Fourier Basis**,
*Journal of Chemical Theory and Computation*.

DOI: https://doi.org/10.1021/acs.jctc.3c01317

Used as context for the exact 2D benchmark architecture:

- local-diabatic treatment of CI wavepacket dynamics;
- Strang splitting;
- Fourier nuclear basis;
- electronic coherence and geometric phase without singular derivative couplings.

## Ben-Nun, Quenneville, and Martínez 2000 — AIMS

M. Ben-Nun, J. Quenneville, and T. J. Martínez,
**Ab Initio Multiple Spawning: Photochemistry from First Principles Quantum Molecular
Dynamics**,
*Journal of Physical Chemistry A* **104**, 5161-5175 (2000).

DOI: https://doi.org/10.1021/jp994174i

Retained as the foundational reference for the AIMS-style coupled Gaussian and
adaptive-spawning framework.

## Curchod, Glover, and Martínez 2020 — SSAIMS

B. F. E. Curchod, W. J. Glover, and T. J. Martínez,
**SSAIMS—Stochastic-Selection Ab Initio Multiple Spawning for Efficient Nonadiabatic
Molecular Dynamics**,
*Journal of Physical Chemistry A* **124**, 6133-6143 (2020).

DOI: https://doi.org/10.1021/acs.jpca.0c04113

Used for modern context on:

- coupled frozen-Gaussian TBFs;
- coherence in AIMS;
- adaptive basis growth and the cost of unnecessary TBFs;
- the formal complete-basis/exact-integral limit.

## Yang, Coe, Kaduk, and Martínez 2009 — optimal spawning

S. Yang, J. D. Coe, B. Kaduk, and T. J. Martínez,
**An "optimal" spawning algorithm for adaptive basis set expansion in nonadiabatic
dynamics**,
*Journal of Chemical Physics* **130**, 134113 (2009).

DOI: https://doi.org/10.1063/1.3103930

Retained for the v0.11/v0.12 spawning design principle:

- spawned functions should add useful coupled basis directions;
- parent-child placement is not arbitrary;
- fewer, better placed spawned functions can be more efficient than uncontrolled
  basis proliferation.

v0.12's projected-bank benchmark additionally demonstrates that crossing a coupling
threshold is not, by itself, evidence that another TBF will reduce the measured
wavefunction or density-matrix residual.


# v0.13 references: adaptive Gaussian basis quality and convergence

## SSAIMS — stochastic selection and adaptive TBF control

B. F. E. Curchod, W. J. Glover, and T. J. Martínez,
**SSAIMS—Stochastic-Selection Ab Initio Multiple Spawning for Efficient Nonadiabatic
Molecular Dynamics**,
*Journal of Physical Chemistry A* **124**, 6133-6143 (2020).

DOI: https://doi.org/10.1021/acs.jpca.0c04113

Used for literature context on:

- AIMS/FMS as an adaptive basis of coupled frozen Gaussian TBFs;
- the formally exact complete-basis/exact-integral limit;
- the computational problem created by uncontrolled TBF growth;
- the need for principled basis selection/removal.

The v0.13 Hilbert-residual and TDSE-defect scores are repository-specific derivations
and are not attributed to SSAIMS.

## Legion — modern modular multiple-spawning implementation

**Legion: A Platform for Gaussian Wavepacket Nonadiabatic Dynamics**,
*Journal of Chemical Theory and Computation*.

DOI: https://doi.org/10.1021/acs.jctc.4c01697

Used for current implementation context on:

- automatic Gaussian basis expansion in strong nonadiabatic regions;
- energy-conserving child creation;
- trajectory elimination/basis control;
- systematic improvement of Gaussian wavepacket calculations by increasing and
  refining the trajectory basis.

The v0.13 residual-driven selection criterion is not claimed to be part of Legion.

## G-MCTDH — variational Gaussian convergence

I. Burghardt, G. A. Worth, and co-workers,
**Multimode quantum dynamics using Gaussian wavepackets: The Gaussian-based
multiconfiguration time-dependent Hartree (G-MCTDH) method applied to the absorption
spectrum of pyrazine**,
*Journal of Chemical Physics* **129**, 174104 (2008).

DOI: https://doi.org/10.1063/1.2996349

Used for contrast and convergence context:

- Gaussian wavepacket bases can be systematically converged toward exact quantum
  dynamics;
- variational/nonclassical motion of the Gaussian basis can be important for
  convergence.

v0.13 does not implement G-MCTDH or vMCG parameter equations; it keeps classically
guided Gaussian motion and improves basis selection through explicit residuals.

## Optimal spawning

S. Yang, J. D. Coe, B. Kaduk, and T. J. Martínez,
**An "optimal" spawning algorithm for adaptive basis set expansion in nonadiabatic
dynamics**,
*Journal of Chemical Physics* **130**, 134113 (2009).

DOI: https://doi.org/10.1063/1.3103930

Retained for the complementary principle that newly created TBFs should be placed so
that they contribute useful parent-child coupling while satisfying energy constraints.

v0.13 differs by ranking missing Hilbert/Galerkin directions directly from a measured
residual.


# v0.14 reference note

v0.14 does not introduce a new literature attribution for its repository-specific
TDSE-defect hysteresis controller, exact leave-one-out pruning formula, complexity
ledger, or Hermitian half-build optimization.

The scientific context remains the AIMS/FMS, SSAIMS, optimal-spawning,
local-diabatic, and variational-Gaussian literature listed above.

Algorithmic-complexity expressions in `docs/releases/v0.14/V14_ALGORITHM_COMPLEXITY.md` are direct
operation-count analyses of the implementation and standard dense linear-algebra/FFT
scaling, not quotations from those papers.


# v0.15 implementation note

v0.15 does not attribute its Gaussian pair-cache, incremental matrix-update, or
cost-aware utility formulas to an external publication.

Those are repository-specific numerical-architecture derivations built on the
Gaussian/AIMS/local-diabatic context already cited in earlier release sections.

The complexity expressions in `docs/releases/v0.15/V15_ALGORITHM_COMPLEXITY.md` are direct operation-count
analyses of the implementation and standard dense-linear-algebra/FFT scaling.


# v0.16 implementation note

The v0.16 KD-tree locality screen, overlap upper-bound derivation, sparse graph
hysteresis, sparse matrix audit, and provider-cost utility are repository-specific
numerical methods.

The release does not attribute those implementation formulas to an external paper.

The broader scientific context remains the Gaussian/AIMS/local-diabatic and
gauge-transport literature already cited in earlier sections.

Complexity statements in `docs/releases/v0.16/V16_ALGORITHM_COMPLEXITY.md` are direct implementation
operation counts and standard sparse/dense linear-algebra considerations.


# v0.17 implementation note

The S/H/T edge score, local omitted-score L2 proxy, one-sided dense-audit relaxation,
and the v0.17 construction-scaling campaign are repository-specific numerical methods.

No external publication is cited as the source of those formulas.

The broader AIMS/FMS, Gaussian, local-diabatic, gauge-transport, and electronic
structure context remains covered by the references listed in earlier release
sections.


# v0.18 implementation note

The v0.18 phase-aligned grid metrics, physical-time adaptive-cadence normalization,
batched residual ranking, sampled omitted-edge audit policy, and multi-axis convergence
campaign are repository-specific numerical methods.

The second-order self-convergence diagnostic is the standard Richardson idea applied
to successive phase-aligned Gaussian wavefunctions.

No external publication is claimed as the source of the repository-specific sampled
audit or batching policy. Earlier Gaussian/AIMS/FMS, gauge-transport, and electronic
structure references remain applicable to the underlying scientific context.


# v0.19 implementation note

The nearest-anchor molecular tracking policy, bounded fallback policy, normalized
provider-cost interface, and synthetic Cartesian LVC molecular backend are
repository-specific implementation choices.

The polynomial state-assignment implementation uses the standard linear-sum assignment
problem for maximum-overlap matching, with additional constrained assignments to
recover the exact second-best score margin.

The PySCF bridge reuses the repository's previously documented SA-CASSCF and
many-electron cross-geometry overlap machinery. No new real PySCF runtime result is
claimed in v0.19.


# v0.20 implementation note

The buffered KD-tree electronic cache, molecular S/H/T edge score, local molecular
omitted-score budget, sampled molecular search-floor controller, and independent
dense-sentinel cache policy are repository-specific numerical methods.

The gauge-covariant pair-centroid transport reuses the discrete electronic gauge-link
formalism documented in earlier releases.

The sparse moving-basis propagation reuses the metric identity

$$
\dot S=T+T^\dagger
$$

and the midpoint/Cayley sparse solve already validated in the analytic sparse releases.

No real PySCF sparse trajectory is claimed in v0.20.


# v0.24.1 static molecular-SOC references

The spin-orbit mean-field context is:

- B. A. Heß, C. M. Marian, U. Wahlgren, and O. Gropen, “A mean-field
  spin-orbit method applicable to correlated wavefunctions,” *Chemical Physics
  Letters* **251**, 365–371 (1996), DOI: 10.1016/0009-2614(96)00119-4.
- F. Neese, “Efficient and accurate approximations to the molecular spin-orbit
  coupling operator and their use in molecular g-tensor calculations,” *Journal of
  Chemical Physics* **122**, 034107 (2005), DOI: 10.1063/1.1829047.
- S. Kotaru, P. Pokhilko, and A. I. Krylov, “Spin-orbit couplings within
  spin-conserving and spin-flipping time-dependent density functional theory:
  implementation and benchmark calculations,” *Journal of Chemical Physics* (2022),
  DOI: 10.1063/5.0130868. This reference documents the same broad combination of a
  Breit-Pauli operator, Wigner--Eckart reduced transition densities, SOMF, and state
  interaction in a different electronic-structure method.

PySCF 2.13.1 supplies the scalar-orbital AO integral, FCI transition-RDM, and
determinant-ladder APIs. The v0.24.1 Clebsch--Gordan finite sum, direct matrix
assembly, capability split, and release gates are repository implementations.
