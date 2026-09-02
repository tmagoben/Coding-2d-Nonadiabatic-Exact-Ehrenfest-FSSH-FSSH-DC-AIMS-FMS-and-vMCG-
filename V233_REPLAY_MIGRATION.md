# v0.23.3 replay format and migration

Replay format 2 uses `molecular_soc_manifest_v233.json` and
`molecular_soc_arrays_v233.npz`. In addition to the component-resolved v0.23.0
data it stores:

- raw directed finite-manifold overlaps;
- independently derived unitary polar transports;
- overlap singular values and the exact quality policy;
- overlap and transport contract identifiers;
- a full derivative-coupling convention and fingerprint;
- a convention-complete provider numerical identity;
- the source legacy fingerprint and, for migrations, an explicit attestation.

Loading recomputes every polar transport from the raw overlap and rejects any
stored mismatch. This prevents a raw contraction from being relabelled as a
transport. Manifest and array digests remain deterministic and corruption
sensitive.

## Migration policy

Format-1 input is never loaded as format 2 implicitly. Migration requires a
typed `LegacyReplayMigrationAttestationV233` bound to the exact legacy dataset.
Accepted dispositions are:

- `not_pyscf_derived`, with the analytic/native internal NAC convention;
- `verified_corrected_v232`, with the certified PySCF 2.13.1 mapping.

`unknown` and `requires_sign_correction` are quarantined. The migrator does not
guess an index order or modify signs. Such data must be regenerated or reviewed
outside the migration path and then attested from independent evidence.

Direct v0.23.3 capture writes format 2 without pretending to be a legacy
migration.
