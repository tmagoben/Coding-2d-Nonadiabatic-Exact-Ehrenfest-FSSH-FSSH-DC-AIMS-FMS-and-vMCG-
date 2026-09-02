import json
from pathlib import Path

from gaussian_dynamics import __version__


def test_v0241_release_metadata_evidence_and_claim_boundary_are_consistent():
    root = Path(__file__).resolve().parents[1]
    pyproject = (root / "pyproject.toml").read_text(encoding="utf-8")
    citation = (root / "CITATION.cff").read_text(encoding="utf-8")
    readme = (root / "README.md").read_text(encoding="utf-8")
    campaign = json.loads(
        (root / "results/v0241_pyscf_static_soc_campaign.json").read_text(
            encoding="utf-8"
        )
    )
    evidence = json.loads(
        (root / "results/v0241_pyscf_static_soc_evidence.json").read_text(
            encoding="utf-8"
        )
    )

    assert __version__ == "0.27.0"
    assert 'version = "0.27.0"' in pyproject
    assert 'pyscf = ["pyscf==2.13.1"]' in pyproject
    assert 'dependencies = ["numpy>=1.24", "scipy>=1.10"]' in pyproject
    assert "sympy" not in pyproject.lower()
    assert "prism" not in pyproject.lower()
    assert "version: 0.27.0" in citation
    assert "date-released: 2026-08-25" in citation
    assert "Current release: v0.27.0" in readme
    assert '[tool.setuptools]\npackages = ["gaussian_dynamics"]' in pyproject

    required = (
        "V241_RELEASE_NOTES.md",
        "V241_PYSCF_STATE_INTERACTION_SOC.md",
        "V241_PROGRAM_ARCHITECTURE.md",
        "V241_ALGORITHM_COMPLEXITY.md",
        "V241_VALIDATION.md",
        "V241_BUILD_VALIDATION.md",
        "docs/18_PYSCF_STATIC_MOLECULAR_SOC.md",
        "examples/127_recompute_v0241_pyscf_static_soc.py",
        "examples/128_recompute_v0241_campaign.py",
        "requirements-pyscf-v241-linux-x86_64-py312.txt",
        "results/v0241_pyscf_static_soc_evidence.json",
        "results/v0241_pyscf_static_soc_campaign.json",
        "gaussian_dynamics/pyscf_state_interaction_soc_v241.py",
        "gaussian_dynamics/pyscf_soc_runtime_v241.py",
        "gaussian_dynamics/v241_benchmark.py",
    )
    assert all((root / name).is_file() for name in required)

    acceptance = campaign["acceptance"]
    claims = campaign["claims"]
    assert acceptance["passed"] is True
    assert acceptance["inherited_gate_count"] == 256
    assert acceptance["runtime_gate_count"] == 39
    assert acceptance["core_gate_count"] == 20
    assert acceptance["new_gate_count"] == 59
    assert acceptance["total_gate_count"] == 315
    assert len(acceptance["checks"]) == 315
    assert all(type(value) is bool for value in acceptance["checks"].values())
    assert all(acceptance["checks"].values())

    assert claims["real_PySCF_BP_SOMF_execution_validated"] is True
    assert claims["direct_molecular_SOC_elements_returned"] is True
    assert claims["static_molecular_SOC_tier_validated"] is True
    assert claims["trajectory_ready_molecular_SOC_validated"] is False
    assert claims["live_molecular_SOC_backend_admitted"] is False
    assert claims["physical_SOC_derivatives_validated"] is False
    assert claims["cross_geometry_SOC_tracking_validated"] is False
    assert claims["ab_initio_SOC_accuracy_validated"] is False
    assert claims["external_molecular_SOC_snapshot_admitted"] is False
    assert claims["native_openmolcas_numeric_crosscheck_implemented"] is False
    assert claims["Prism_runtime_dependency_required"] is False

    assert evidence["schema"] == "gnd-pyscf-static-soc-runtime-evidence-v0.24.1"
    assert evidence["audit"]["passed"] is True
    assert len(evidence["audit"]["checks"]) == 39
    assert evidence["result"]["capabilities"]["tier"] == "static_soc"
    assert evidence["result"]["trajectory_ready"] is False
    assert evidence["result"]["matrices"]["state_order"] == [
        "D1(M=+1/2)",
        "D1(M=-1/2)",
        "D2(M=+1/2)",
        "D2(M=-1/2)",
        "D3(M=+1/2)",
        "D3(M=-1/2)",
    ]
    assert "H_soc" in evidence["result"]["matrices"]
    assert evidence["audit"]["metrics"]["H_soc_frobenius_norm_cm_inverse"] > 100
