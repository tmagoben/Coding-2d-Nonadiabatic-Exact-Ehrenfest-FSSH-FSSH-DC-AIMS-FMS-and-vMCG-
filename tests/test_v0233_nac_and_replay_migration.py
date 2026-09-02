from dataclasses import replace
import json

import numpy as np
import pytest

from gaussian_dynamics.molecular_soc_replay_v233 import (
    FileBackedMolecularSOCProviderV233,
    load_molecular_soc_replay_v233,
    migrate_molecular_soc_replay_v230_to_v233,
)
from gaussian_dynamics.nac_compatibility_v233 import (
    LegacyReplayMigrationAttestationV233,
    analytic_nac_convention_v233,
    corrected_pyscf_nac_convention_v233,
    require_snapshot_nac_identity_v233,
)
from gaussian_dynamics.v230_benchmark import (
    build_v230_reference_replay,
    v230_reference_coordinates,
)


def _analytic_attestation(dataset, disposition="not_pyscf_derived"):
    return LegacyReplayMigrationAttestationV233(
        legacy_dataset_fingerprint=dataset.dataset_fingerprint,
        nac_disposition=disposition,
        evidence="deterministic analytic provider emits internal d[i,j] directly",
    )


def test_corrected_pyscf_full_overlap_and_etf_identities_are_distinct():
    full = corrected_pyscf_nac_convention_v233()
    etf = corrected_pyscf_nac_convention_v233(use_etfs=True)
    assert full.fingerprint() != etf.fingerprint()
    assert not full.use_etfs
    assert etf.use_etfs

    assert require_snapshot_nac_identity_v233(
        {"v233_nac_convention_fingerprint": full.fingerprint()}, full
    )
    with pytest.raises(ValueError, match="lacks v0.23.3"):
        require_snapshot_nac_identity_v233({}, full)
    with pytest.raises(ValueError, match="identity mismatch"):
        require_snapshot_nac_identity_v233(
            {"v233_nac_convention_fingerprint": etf.fingerprint()}, full
        )


@pytest.mark.parametrize("disposition", ["unknown", "requires_sign_correction"])
def test_legacy_unknown_or_wrong_sign_data_are_quarantined(tmp_path, disposition):
    legacy = build_v230_reference_replay(tmp_path / "legacy")
    attestation = _analytic_attestation(legacy, disposition=disposition)
    with pytest.raises(ValueError, match="quarantined"):
        migrate_molecular_soc_replay_v230_to_v233(
            legacy.manifest_path,
            tmp_path / "v233",
            nac_convention=analytic_nac_convention_v233(),
            migration_attestation=attestation,
        )


def test_attested_replay_migration_is_deterministic_and_versioned(tmp_path):
    legacy = build_v230_reference_replay(tmp_path / "legacy")
    convention = analytic_nac_convention_v233()
    attestation = _analytic_attestation(legacy)
    first = migrate_molecular_soc_replay_v230_to_v233(
        legacy.manifest_path,
        tmp_path / "first",
        nac_convention=convention,
        migration_attestation=attestation,
    )
    second = migrate_molecular_soc_replay_v230_to_v233(
        legacy.manifest_path,
        tmp_path / "second",
        nac_convention=convention,
        migration_attestation=attestation,
    )

    assert first.dataset_fingerprint == second.dataset_fingerprint
    assert first.manifest_path.read_bytes() == second.manifest_path.read_bytes()
    assert first.arrays_path.read_bytes() == second.arrays_path.read_bytes()
    manifest = json.loads(first.manifest_path.read_text(encoding="utf-8"))
    assert manifest["format_version"] == 2
    assert manifest["release"] == "v0.23.3"
    assert manifest["nac_convention_fingerprint"] == convention.fingerprint()
    assert "overlap_quality_policy" in manifest
    assert "transport_contract" in manifest


def test_v233_loader_rejects_unmigrated_v230_manifest(tmp_path):
    legacy = build_v230_reference_replay(tmp_path / "legacy")
    with pytest.raises(ValueError, match="quarantined"):
        load_molecular_soc_replay_v233(legacy.manifest_path)


def test_file_provider_exposes_raw_overlap_and_transport_separately(tmp_path):
    legacy = build_v230_reference_replay(tmp_path / "legacy")
    convention = analytic_nac_convention_v233()
    migrated = migrate_molecular_soc_replay_v230_to_v233(
        legacy.manifest_path,
        tmp_path / "v233",
        nac_convention=convention,
        migration_attestation=_analytic_attestation(legacy),
    )
    provider = FileBackedMolecularSOCProviderV233(migrated.manifest_path)
    coordinates = v230_reference_coordinates()
    left = provider.evaluate_snapshot(coordinates[0])
    right = provider.evaluate_snapshot(coordinates[1])

    raw = provider.snapshot_overlap(left, right)
    transport = provider.snapshot_transport(left, right)
    assert np.allclose(transport.conj().T @ transport, np.eye(4))
    assert left.metadata["v233_nac_convention_fingerprint"] == (
        convention.fingerprint()
    )
    assert np.allclose(raw, migrated.overlaps[0, 1])
    assert np.allclose(transport, migrated.overlap_transports[0, 1])


def test_stored_raw_overlap_cannot_be_substituted_for_transport(tmp_path):
    legacy = build_v230_reference_replay(tmp_path / "legacy")
    migrated = migrate_molecular_soc_replay_v230_to_v233(
        legacy.manifest_path,
        tmp_path / "v233",
        nac_convention=analytic_nac_convention_v233(),
        migration_attestation=_analytic_attestation(legacy),
    )
    bad = np.asarray(migrated.overlap_transports, dtype=complex).copy()
    bad[0, 1] = 0.9 * np.eye(bad.shape[-1])
    corrupted = replace(migrated, overlap_transports=bad)
    with pytest.raises(ValueError, match="transport differs"):
        corrupted.validate()
