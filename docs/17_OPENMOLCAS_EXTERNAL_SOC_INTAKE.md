# OpenMolcas external-SOC intake

v0.24.0 adds a method-specific evidence plane above the representation-neutral
electronic operator contract. Its core rule is simple: numerical plausibility is not
provenance. A matrix can be Hermitian, time-reversal symmetric, converged under finite
differences, and still not be an ab-initio result.

The intake therefore combines numerical gates with native artifact signatures,
content digests, exact source classification, an out-of-band trust policy, and an
independent accuracy dossier. The bundled generator creates only a protocol fixture;
its non-execution marker and non-HDF5 placeholder make its status machine-checkable.

See `V240_OPENMOLCAS_PROTOCOL.md`, `V240_EXTERNAL_SNAPSHOT_ADMISSION.md`, and
`V240_PROGRAM_ARCHITECTURE.md` for the complete protocol, gate sequence, and diagram.

