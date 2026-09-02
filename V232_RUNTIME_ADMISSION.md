# v0.23.2 trust-anchored runtime admission

v0.23.2 closes several ways that plausible metadata could previously cross the
future molecular-SOC admission boundary without proving a real implementation.

## Admission order

1. Validate object structure and require callable engine methods.
2. Compare the complete backend/method identity with a caller-supplied trusted
   identity: backend/version, source kind, adapter, method, basis, active space,
   SOC operator, scalar-relativistic treatment, derivative method, and NAC
   convention.
3. Require all declared capabilities.
4. Require exactly seven convergence stages: SCF, correlated wavefunction,
   state-interaction SOC, spin-free gradients, SOC derivatives, derivative
   connections, and many-electron overlaps.
5. Bind execution challenge, runtime probe, environment, raw inventory, receipts,
   dossier, replay, and parsed results.
6. Execute the exact trusted validator type and demand typed parser/execution proof.
7. Run the inherited molecular physical, derivative, and symmetry audits.

Point and snapshot metadata cannot be merged to manufacture convergence. Legacy
five-stage vocabularies, duplicate convergence namespaces, non-callable method
placeholders, self-asserted validator identities, unattested live sources, and
Boolean-only parser claims are rejected.

An external snapshot may be parsed without pretending that it was freshly
executed. A live admission additionally requires a fresh execution bound to the
challenge and runtime evidence. The two paths are intentionally distinct.

The negative-control campaign proves that these barriers fail closed. No external
or live molecular-SOC source is admitted in v0.23.2.
