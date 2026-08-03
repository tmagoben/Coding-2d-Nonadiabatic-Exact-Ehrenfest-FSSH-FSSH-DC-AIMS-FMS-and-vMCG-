# Validation and Benchmark Results

All automated tests and all sequential examples were executed successfully.

## Static 1D electronic problem

- Maximum difference between Hellmann--Feynman and finite-difference
  derivative couplings: \(3.693182\times10^{-5}\).
- Maximum error after reconstructing the original diabatic potential:
  \(1.605347\times10^{-8}\).

## 1D dynamics

The \(256\)-point periodic-grid calculation was propagated to
\(t=1000\) atomic units.

- Maximum FFT norm error: \(7.194245\times10^{-14}\).
- Maximum direct-propagation norm error: \(1.776357\times10^{-15}\).
- Final FFT/direct fidelity: numerically \(1.000000000000\).
- Final phase-aligned wavefunction error:
  \(3.012624\times10^{-7}\).
- Final FFT adiabatic populations:
  \[
  (P_0,P_1)=(0.39488324,0.60511676).
  \]
- Final direct adiabatic populations:
  \[
  (P_0,P_1)=(0.39488328,0.60511672).
  \]

## Static 2D conical intersection

- Maximum analytic versus Hellmann--Feynman derivative-coupling error:
  \(2.664535\times10^{-15}\).
- Closed-loop derivative-coupling integral:
  \[
  3.141592653590\approx\pi.
  \]
- ADT holonomy:
  \[
  A_fA_i^\mathsf{T}\approx-I.
  \]

## 2D dynamics

- Maximum norm error in the \(64\times64\) FFT calculation:
  \(3.841372\times10^{-13}\).
- Final diabatic populations:
  \[
  (0.40881917,0.59118083).
  \]
- Final adiabatic populations:
  \[
  (0.71136889,0.28863111).
  \]
- Tiny-grid FFT/direct benchmark fidelity:
  \(1.000000000000\).
- Tiny-grid phase-aligned error:
  \(6.477682\times10^{-9}\).

## Strang time-step convergence

The errors for

\[
\Delta t=(0.8,0.4,0.2,0.1)
\]

were

\[
(5.271974\times10^{-8},
 1.317981\times10^{-8},
 3.294945\times10^{-9},
 8.237350\times10^{-10}).
\]

A log--log fit gives the global order

\[
p=2.000006,
\]

confirming the expected \(O(\Delta t^2)\) global error of the symmetric
split-operator method.
