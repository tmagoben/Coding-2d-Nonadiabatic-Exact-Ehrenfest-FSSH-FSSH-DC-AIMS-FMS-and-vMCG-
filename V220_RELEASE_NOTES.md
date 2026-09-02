# v0.22.0 release notes

v0.22.0 is the first release in this project to introduce a nonzero physical
spin-orbit Hamiltonian and its physical nuclear derivative. The implementation is
deliberately limited to exactly reproducible analytic models so that the operator
contract can be validated before any electronic-structure backend is trusted with
SOC data.

## Added

- An even-electron model containing one singlet and one complete triplet in the
  ordered basis
  \(\{|S\rangle,|T_{-1}\rangle,|T_0\rangle,|T_{+1}\rangle\}\).
- An odd-electron model containing two complete Kramers doublets in the ordered basis
  \(\{|D_1,+\tfrac12\rangle,|D_1,-\tfrac12\rangle,
  |D_2,+\tfrac12\rangle,|D_2,-\tfrac12\rangle\}\).
- Explicit decompositions

  \[
  H=H_{\rm sf}+H_{\rm SOC},\qquad
  K_a=K_{a,{\rm sf}}+K_{a,{\rm SOC}},
  \]

  with analytic physical derivatives in atomic units.
- Even- and odd-electron time-reversal conventions, including the correct
  transformation of the antiunitary representation under a general complex gauge.
- Gauge-covariant singlet/triplet and doublet-root projectors.
- An independent one-dimensional vector split-operator grid propagator.
- SOC-active dense, sparse, moving-frame, and deterministic-restart validation.
- A canonical 53-gate campaign: all 21 v0.21.4 gates plus 32 new physical-SOC gates.

## Acceptance result

The canonical campaign passes **53/53** gates. It covers zero, constant, and
coordinate-dependent SOC; physical-force finite differences; time reversal; Kramers
degeneracy; arbitrary complex gauges; exact-grid norm, energy, resolution, box, and
timestep convergence; short-time Gaussian/grid population agreement; deterministic
restart; provenance mismatch; and checkpoint corruption.

The cumulative source regression suite passes **302 tests**.

## Compatibility and boundaries

- Spin-free dynamics remains a permanent first-class mode.
- Even- and odd-electron model spaces are separate and are never mixed.
- The v0.21.4 archive is not modified by this release.
- No ab-initio SOC backend, PySCF SOC trajectory, molecular SOC accuracy, external
  magnetic field, or production AIMS claim is made.

The detailed equations are in `V220_DERIVATIONS.md`; measured validation results are
in `V220_VALIDATION.md`; the stable public contract and claim boundary are summarized
in `docs/11_PHYSICAL_ANALYTIC_SOC.md`.
