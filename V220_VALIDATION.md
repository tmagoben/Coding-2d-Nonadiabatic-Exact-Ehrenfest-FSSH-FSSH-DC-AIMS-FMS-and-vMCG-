# v0.22.0 validation record

## Acceptance structure

The canonical release campaign contains **53 gates**:

- 21 inherited v0.21.4 gates, unchanged;
- 8 shared physical-SOC contract gates;
- 6 singlet–triplet gates;
- 8 doublet/Kramers gates;
- 5 independent exact-grid and convergence gates;
- 5 SOC restart and failure-control gates.

All **53/53** gates pass in the release environment.

The cumulative source suite, including all inherited releases, reports:

```text
302 passed in 124.68s
```

## Physical operator and symmetry checks

| Quantity | Singlet–triplet | Two doublets |
|---|---:|---:|
| Cross-geometry \(K\) differential error | 1.8665508646787508e-15 | 2.8189256543263856e-15 |
| SOC-force finite-difference error | 3.108180623686077e-15 | within the 2.0e-10 gate |
| Maximum Kramers pair splitting | not applicable | 4.336808689942018e-18 |

The campaign additionally requires exact zero-SOC operator equivalence, a nonzero
constant SOC Hamiltonian with exactly zero SOC derivative, correct even- and
odd-electron time-reversal squares, complete projectors, and rejection of a Hermitian
but physically wrong SOC derivative.

## Independent exact-grid propagation

For \(\Delta t=0.04\), 100 steps, and a 256-point periodic grid on \([-8,8)\):

| Model | Maximum norm drift | Maximum energy drift | Final transferred population |
|---|---:|---:|---:|
| Singlet–triplet | 4.440892098500626e-15 | 6.350649173203493e-14 | 1.883844585186722e-04 triplet |
| Two doublets | 1.021405182655144e-14 | 1.1173440644940413e-13 | 7.362475396071825e-05 in doublet 2 |

The timestep differences are

\[
\epsilon_{0.08,0.04}=1.5208404669523864\times10^{-10},
\qquad
\epsilon_{0.04,0.02}=3.802085470115884\times10^{-11},
\]

giving observed order

\[
p=2.0000059562889683.
\]

The maximum transferred-population change under grid-spacing and box refinement is
\(8.632959798415829\times10^{-18}\).

## Gaussian versus grid populations

At the deliberately short comparison time \(t=0.2\), where the single-Gaussian
representation error remains controlled:

| Model | Gaussian population | Grid population | Absolute error |
|---|---:|---:|---:|
| Singlet–triplet | 4.636143889725949e-07 | 4.710200388641343e-07 | 7.405649891539424e-09 |
| Two doublets | 1.819251440398155e-07 | 1.8417504618510336e-07 | 2.2499021452878545e-09 |

Both satisfy the release threshold of \(10^{-8}\). This gate compares physical
projector populations, not gauge-dependent amplitudes.

## Gauge and restart checks

- Moving-complex-frame Gaussian dynamics coefficient error:
  \(5.936673241766072\times10^{-16}\); position and momentum errors are zero.
- Dense SOC restart coefficient error:
  \(1.3877787809430044\times10^{-17}\); position and momentum errors are zero.
- Sparse SOC restart coefficient error:
  \(1.6531559496939515\times10^{-16}\); position and momentum errors are zero.
- Sparse edge \([[3,8]]\) is preserved across the restart boundary.
- Changed SOC provenance and deliberately corrupted checkpoint data are rejected.

## Failure controls and claim boundary

Passing controls include an incomplete-doublet rejection, an electron-parity check, a
broken-Kramers fixture, a wrong-but-Hermitian SOC derivative, an untransformed
time-reversal representation in a complex gauge, and changed restart provenance.

These results validate the implemented analytic models and the framework contracts.
They do not validate molecular SOC matrix elements, an ab-initio derivative method, a
PySCF SOC runtime, external magnetic fields, or molecular predictive accuracy.
