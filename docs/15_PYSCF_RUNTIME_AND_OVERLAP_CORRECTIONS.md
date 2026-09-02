# PySCF runtime and overlap corrections

v0.23.2 turns the optional PySCF boundary into a validated real spin-free runtime
without overstating molecular-SOC support.

## Three corrected contracts

1. **NAC orientation.** Internal `d[i,j]=<Phi_i|d Phi_j>` uses PySCF 2.13.1
   `state=(i,j)` in the empirically certified production adapter. The mapping is
   verified by phase-aligned many-electron overlap central differences.
2. **ETF semantics.** `use_etfs=False` corresponds to the full overlap derivative;
   the ETF quantity is kept separate because it removes translation.
3. **Finite selected manifolds.** Cross-geometry overlaps are contractions. The
   framework requires self identity, reciprocity, and singular values no greater
   than one, rather than exact cross-geometry isometry.

## Runtime evidence

The H3+ SA-CASSCF fixture proves real energies, analytic gradients, NACs, and
cross-geometry many-electron overlaps in a hash-locked PySCF 2.13.1 runtime. The
runtime fingerprint verifies installed package content, not only version strings.

## Admission remains closed

PySCF core does not itself supply the complete method-specific state-interaction
SOC engine required here. A future admitted source must still provide physical
SOC derivatives, complete multiplets, exact method identity, seven convergence
stages, raw receipts and artifacts, replay/dossier binding, and trusted typed
parser/execution proof. v0.23.2 strengthens and tests that boundary but admits no
external or live molecular-SOC provider.

See the root-level v0.23.2 release, runtime, erratum, overlap, admission,
validation, complexity, and architecture documents for the complete record.
