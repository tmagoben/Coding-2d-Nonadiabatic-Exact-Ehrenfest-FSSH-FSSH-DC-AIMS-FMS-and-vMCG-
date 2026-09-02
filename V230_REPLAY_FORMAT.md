# v0.23.0 molecular SOC replay format

## Files

A replay dataset is one directory containing exactly the named format files:

```text
molecular_soc_manifest_v230.json
molecular_soc_arrays_v230.npz
```

The canonical JSON manifest uses sorted keys and compact separators. The NPZ writer
sorts array names and fixes every ZIP timestamp to 1980-01-01, so the same logical data
produce identical bytes.

## Stored arrays

For $R$ records, $d$ nuclear coordinates, and $s$ electronic states, the archive stores:

| Array | Shape | Meaning |
|---|---:|---|
| `q` | `(R,d)` | exact generalized coordinates in bohr |
| `H_spin_free` | `(R,s,s)` | spin-free Hamiltonians in hartree |
| `K_spin_free` | `(R,d,s,s)` | physical spin-free derivatives |
| `H_soc` | `(R,s,s)` | SOC Hamiltonians in hartree |
| `K_soc` | `(R,d,s,s)` | physical SOC derivatives |
| `connection_q` | `(R,d,s,s)` | derivative connections |
| `mass_matrix_q_au` | `(R,d,d)` | nuclear mass matrices in atomic units |
| `overlaps` | `(R,R,s,s)` | every ordered cross-record electronic overlap |
| `converged` | `(R,)` | per-record electronic convergence flags |
| `time_reversal_matrix` | `(s,s)` | numerical antiunitary representation |
| `projectors` | `(P,s,s)` | named physical projectors |

The manifest freezes dimensions, state-sector metadata, projector names, the complete
molecular admission contract, inherited operator provenance, their fingerprints, the
source-provider fingerprint, the array SHA-256, and a dataset fingerprint.

## Integrity rules

Loading rejects:

- an unknown format or version;
- an unexpected array member set, shape, dtype role, or non-finite value;
- an array SHA-256 mismatch;
- a manifest, provenance, or contract fingerprint mismatch;
- coordinate collisions under the declared digit policy;
- non-Boolean convergence flags or contract/record convergence disagreement;
- overlap identity, isometry, or reciprocity failure;
- inconsistent parity, time reversal, projector, or dimension data.

Replay tokens carry the dataset fingerprint. Cross-dataset overlap requests are
rejected even when record indices happen to match.

## Exact-record policy

The file-backed provider accepts only coordinates explicitly stored in `q`, rounded
under the manifest's coordinate-digit policy. Interpolation and extrapolation are
forbidden because they would invent electronic operators, derivatives, and wavefunction
overlaps that were never calculated.

## Reference datasets

v0.23.0 includes two nine-record, one-coordinate, four-state fixtures:

- even singlet–triplet fingerprint:
  `a0c90420ace96c899b5e033a8b03d43cd316eb17511985140541be4e4dec8255`;
- odd two-Kramers-doublet fingerprint:
  `9e12baf2950ee2f7fe6ec5debe183cc3c5b2ff1e25ef3df20fe636e05a2c9acb`.

They validate the format and both electron-parity paths. Their source kind is
`validation_fixture`; they are not molecular accuracy evidence.
