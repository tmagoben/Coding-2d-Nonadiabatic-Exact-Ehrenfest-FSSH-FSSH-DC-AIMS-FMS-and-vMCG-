# Corrective SOC hardening

v0.22.1 is the admission barrier between the analytic SOC milestone and a future
molecular backend. It adds no new physical Hamiltonian; it strengthens how an SOC
provider is trusted.

The central changes are:

1. certify transported full-matrix derivatives of `H_spin_free` and `H_soc` separately
   for every nuclear coordinate;
2. admit only a single electron-parity and charge sector with complete multiplets;
3. require both the time-reversal square and time-reversal unitarity;
4. bind numerical time-reversal and projector data into provenance identity;
5. keep the exact-grid oracle fixed-frame and obtain mass from the emitted operator
   contract;
6. require the exact final sample and precompute static split operators;
7. add physical SOC convergence with respect to Gaussian basis size and sparse graph
   threshold.

The adversarial derivative fixture is essential. If

\[
K_{\mathrm{sf}}\mapsto K_{\mathrm{sf}}-\Delta,\qquad
K_{\mathrm{SOC}}\mapsto K_{\mathrm{SOC}}+\Delta,
\]

then total K remains unchanged. If the sampled density also satisfies
`Tr(rho Delta)=0`, a scalar-force check remains unchanged. Separate matrix-valued
finite differences are therefore required to expose both errors.

Likewise, testing only `JJ* = +/-I` does not prove that the antiunitary representation
is norm preserving. The independent condition `J†J = I` is mandatory.

The full equations, measured residuals, and cost model are in
`V221_CORRECTIVE_HARDENING.md`, `V221_VALIDATION.md`, and
`V221_ALGORITHM_COMPLEXITY.md`.

No molecular SOC backend is admitted in v0.22.1. A v0.23 backend must satisfy these
contracts and add independent molecular reference and method/basis convergence evidence.
