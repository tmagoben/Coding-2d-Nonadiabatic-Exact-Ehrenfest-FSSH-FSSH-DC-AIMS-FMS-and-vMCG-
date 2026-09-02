from dataclasses import replace
from pathlib import Path

import pytest

from gaussian_dynamics import (
    RawArtifactRecordV231,
    audit_molecular_soc_provider_v230,
    audit_molecular_soc_provider_v231,
    build_v231_admission_bundle,
    load_molecular_soc_dossier_v231,
    probe_pyscf_soc_adapter_v231,
    require_pyscf_soc_adapter_v231,
)


def test_even_fixture_passes_raw_protocol_but_not_real_admission(tmp_path):
    bundle = build_v231_admission_bundle(tmp_path / "even")
    protocol = audit_molecular_soc_provider_v231(
        bundle["provider"], bundle["center"], bundle["dossier_path"],
        requirement="protocol",
    )
    real = audit_molecular_soc_provider_v231(
        bundle["provider"], bundle["center"], bundle["dossier_path"],
        requirement="real",
    )

    assert protocol.protocol_passed
    assert protocol.passed
    assert not real.real_backend_admitted
    assert not real.passed
    assert real.source_kind == "validation_fixture"


def test_odd_doublet_fixture_uses_connected_manifold_evidence(tmp_path):
    bundle = build_v231_admission_bundle(tmp_path / "odd", odd=True)
    report = audit_molecular_soc_provider_v231(
        bundle["provider"], bundle["center"], bundle["dossier_path"],
        requirement="protocol",
    )

    assert report.protocol_passed
    assert report.tracking_report["minimum_overlap"] == 1.0
    assert report.tracking_report["minimum_margin"] == 1.0
    assert bundle["provider"].soc_symmetry_contract.electron_parity == "odd"


def test_dossier_roundtrip_and_raw_inventory_are_deterministic(tmp_path):
    first = build_v231_admission_bundle(tmp_path / "first")
    second = build_v231_admission_bundle(tmp_path / "second")

    assert first["dossier_path"].read_bytes() == second["dossier_path"].read_bytes()
    assert first["dossier"].fingerprint() == second["dossier"].fingerprint()
    assert {x.name: x.sha256 for x in first["dossier"].artifacts} == {
        x.name: x.sha256 for x in second["dossier"].artifacts
    }


def test_raw_artifact_corruption_is_rejected(tmp_path):
    bundle = build_v231_admission_bundle(tmp_path / "bundle")
    artifact = bundle["dossier"].artifacts[-1]
    path = bundle["directory"] / artifact.relative_path
    path.write_bytes(path.read_bytes() + b"changed")

    with pytest.raises(ValueError, match="mismatch"):
        load_molecular_soc_dossier_v231(bundle["dossier_path"])


def test_receipt_coordinate_and_convergence_are_independent_gates(tmp_path):
    bundle = build_v231_admission_bundle(tmp_path / "bundle")
    dossier = bundle["dossier"]
    receipts = list(dossier.receipts)
    receipts[0] = replace(receipts[0], q_bohr=(receipts[0].q_bohr[0] + 0.01,))
    with pytest.raises(ValueError, match="coordinate"):
        replace(dossier, receipts=tuple(receipts)).validate(
            bundle_directory=bundle["directory"], dataset=bundle["dataset"]
        )

    receipts[0] = replace(dossier.receipts[0], soc_converged=False)
    report_dossier = replace(dossier, receipts=tuple(receipts))
    report = audit_molecular_soc_provider_v231(
        bundle["provider"], bundle["center"], report_dossier,
        requirement="protocol", bundle_directory=bundle["directory"],
    )
    assert not report.protocol_passed
    assert not report.dossier_protocol_checks["trajectory_receipt_convergence"]


def test_artifact_paths_cannot_escape_bundle():
    record = RawArtifactRecordV231(
        name="escape",
        relative_path="../escape.json",
        role="calculation_output",
        sha256="0" * 64,
        size_bytes=0,
    )
    with pytest.raises(ValueError, match="inside"):
        record.validate()


def test_synthetic_external_relabel_is_blocked_without_parser_attestation(tmp_path):
    bundle = build_v231_admission_bundle(
        tmp_path / "external", source_kind="external_ab_initio_snapshot"
    )
    inherited = audit_molecular_soc_provider_v230(
        bundle["provider"], bundle["center"], require_real_backend=True
    )
    hardened = audit_molecular_soc_provider_v231(
        bundle["provider"], bundle["center"], bundle["dossier_path"],
        requirement="external",
    )

    assert inherited.real_backend_admitted
    assert hardened.protocol_passed
    assert not hardened.external_snapshot_admitted
    assert not hardened.external_admission_checks[
        "backend_artifact_parser_validated"
    ]


def test_pyscf_adapter_probe_fails_closed_without_complete_runtime():
    probe = probe_pyscf_soc_adapter_v231()

    if probe.installed:
        assert not probe.live_admission_ready
        with pytest.raises(RuntimeError):
            require_pyscf_soc_adapter_v231()
    else:
        with pytest.raises(ImportError):
            require_pyscf_soc_adapter_v231()
