# v0.24.2 algorithmic complexity

Let `n` be the AO count, `m` the scalar-MO count, `R` the number of spin-free roots,
`s` the number of complete spin microstates, `N_det` the active-CI determinant count,
and `L` the number of finite-difference step sizes.

## Direct-JK SOMF

The v0.24.1 explicit two-electron route materializes three `n x n x n x n` tensors,
requiring `O(n^4)` storage (equivalently a rank-five array including the three
Cartesian components). v0.24.2 submits three contractions to PySCF's shell-blocked JK
driver and retains only `O(n^2)` AO densities/output matrices per component. The
formal integral work remains implementation- and screening-dependent and can approach
`O(n^4)`, but the production memory bottleneck is removed.

The MO transformation costs `O(3 n^2 m + 3 n m^2)` with standard staged matrix
multiplication, while the state-interaction contraction depends on active-orbital and
root counts rather than the full AO four-index tensor.

## Wavefunction overlap and transport

Restricted-CASSCF root overlaps sum determinant-pair contributions. In the general
dense form this is expensive—up to `O(R^2 N_det^2 n_occ^3)` for determinant-overlap
evaluation—although spin/determinant structure and the small validation active space
reduce the actual work. This is the present scaling target for future optimization.

Lifting root overlaps to a complete microstate matrix is `O(s^2)` storage/work.
The SVD/polar factor, operator transport, and dense residual checks are each
`O(s^3)` time and `O(s^2)` storage.

## Finite-difference scan

With one reusable center and `L` centered step sizes, the scan performs `1+2L`
independent electronic-structure calculations. The canonical v0.24.2 case has
`L=3`, hence seven snapshots. It retains six endpoint compact receipts and three
derivative records. A full Cartesian extension for `N_atom` atoms would require
`1+6N_atom L` independent displaced calculations before symmetry, batching, or
analytic derivatives; v0.24.2 deliberately does not claim that extension.

## Practical next bottlenecks

The production rank-five AO memory problem is solved. The next important scaling
work is therefore:

1. reuse/warm-start orbitals and CI roots across connected geometries;
2. avoid general determinant-pair overlap work where transition-density or sparse-CI
   structure can be exploited;
3. cache reciprocal overlap/polar blocks with provenance-safe keys;
4. introduce analytic or response-theory component derivatives before attempting a
   large full-Cartesian scan.
