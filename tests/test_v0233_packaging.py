import json
from pathlib import Path

from gaussian_dynamics import __version__


def test_inherited_v0240_records_and_current_release_metadata_are_consistent():
    root = Path(__file__).resolve().parents[1]
    pyproject = (root / "pyproject.toml").read_text(encoding="utf-8")
    citation = (root / "CITATION.cff").read_text(encoding="utf-8")
    readme = (root / "README.md").read_text(encoding="utf-8")
    campaign = json.loads(
        (root / "results/v0240_external_soc_intake_campaign.json").read_text(
            encoding="utf-8"
        )
    )

    assert __version__ == "0.27.0"
    assert 'version = "0.27.0"' in pyproject
    assert 'pyscf = ["pyscf==2.13.1"]' in pyproject
    assert "version: 0.27.0" in citation
    assert "date-released: 2026-08-25" in citation
    assert "Current release: v0.27.0" in readme
    assert '[tool.setuptools]\npackages = ["gaussian_dynamics"]' in pyproject

    required = (
        "V240_RELEASE_NOTES.md",
        "V240_OPENMOLCAS_PROTOCOL.md",
        "V240_EXTERNAL_SNAPSHOT_ADMISSION.md",
        "V240_PROGRAM_ARCHITECTURE.md",
        "V240_ALGORITHM_COMPLEXITY.md",
        "V240_BUILD_VALIDATION.md",
        "V240_VALIDATION.md",
        "docs/17_OPENMOLCAS_EXTERNAL_SOC_INTAKE.md",
        "examples/125_create_v0240_protocol_fixture.py",
        "examples/126_recompute_v0240_campaign.py",
        "results/v0240_external_soc_intake_campaign.json",
        "V233_RELEASE_NOTES.md",
        "V233_FINITE_MANIFOLD_TRANSPORT.md",
        "V233_REPLAY_MIGRATION.md",
        "V233_NAC_COMPATIBILITY.md",
        "V233_SOC_CONVENTION.md",
        "V233_RUNTIME_PROFILES.md",
        "V233_VALIDATION.md",
        "V233_ALGORITHM_COMPLEXITY.md",
        "V233_PROGRAM_ARCHITECTURE.md",
        "V233_BUILD_VALIDATION.md",
        "docs/16_TRANSPORT_AND_COMPATIBILITY_HARDENING.md",
        "examples/124_recompute_v0233_campaign.py",
        "requirements-pyscf-v233-linux-x86_64-py312.txt",
        "results/v0233_transport_compatibility_campaign.json",
    )
    assert all((root / name).is_file() for name in required)

    acceptance = campaign["acceptance"]
    claims = campaign["claims"]
    assert acceptance["passed"] is True
    assert acceptance["inherited_gate_count"] == 208
    assert acceptance["new_gate_count"] == 48
    assert acceptance["total_gate_count"] == 256
    assert len(acceptance["checks"]) == 256
    assert all(type(value) is bool for value in acceptance["checks"].values())
    assert all(acceptance["checks"].values())

    assert claims["openmolcas_rassi_so_protocol_frozen"] is True
    assert claims["strict_bundle_artifact_parser_validated"] is True
    assert claims["transported_cartesian_soc_derivative_protocol_validated"] is True
    assert claims["independent_accuracy_evidence_schema_validated"] is True
    assert claims["admission_bound_frozen_snapshot_dynamics_validated"] is True
    assert claims["protocol_fixture_validated"] is True
    assert claims["external_molecular_SOC_snapshot_admitted"] is False
    assert claims["live_molecular_SOC_backend_admitted"] is False
    assert claims["ab_initio_SOC_validated"] is False
    assert claims["openmolcas_runtime_executed"] is False
    assert claims["native_openmolcas_numeric_crosscheck_implemented"] is False

    admission = campaign["admission_controls"]
    assert admission["protocol_fixture_passes_protocol_audit"] is True
    assert admission["protocol_fixture_not_external_admitted"] is True
    assert admission["synthetic_relabel_rejected"] is True
