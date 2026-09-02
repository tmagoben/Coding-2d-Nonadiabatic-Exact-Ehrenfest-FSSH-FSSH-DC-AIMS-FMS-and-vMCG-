import json
from pathlib import Path

from gaussian_dynamics import __version__


def test_v0250_metadata_evidence_and_claim_boundary_are_consistent():
    root = Path(__file__).resolve().parents[1]
    pyproject = (root / "pyproject.toml").read_text(encoding="utf-8")
    citation = (root / "CITATION.cff").read_text(encoding="utf-8")
    readme = (root / "README.md").read_text(encoding="utf-8")
    campaign = json.loads(
        (root / "results/v0250_variational_soc_campaign.json").read_text(
            encoding="utf-8"
        )
    )
    evidence = json.loads(
        (root / "results/v0250_variational_soc_evidence.json").read_text(
            encoding="utf-8"
        )
    )

    assert __version__ == "0.27.0"
    assert 'version = "0.27.0"' in pyproject
    assert 'pyscf = ["pyscf==2.13.1"]' in pyproject
    assert 'dependencies = ["numpy>=1.24", "scipy>=1.10"]' in pyproject
    assert "version: 0.27.0" in citation
    assert "date-released: 2026-08-25" in citation
    assert "Current release: v0.27.0" in readme
    assert '[tool.setuptools]\npackages = ["gaussian_dynamics"]' in pyproject

    required = (
        "V250_RELEASE_NOTES.md",
        "V250_VARIATIONAL_SOC_DYNAMICS.md",
        "V250_INTEGRATOR_DECISION.md",
        "V250_PROGRAM_ARCHITECTURE.md",
        "V250_ALGORITHM_COMPLEXITY.md",
        "V250_VALIDATION.md",
        "V250_BUILD_VALIDATION.md",
        "docs/20_VARIATIONAL_SOC_DYNAMICS.md",
        "examples/131_recompute_v0250_variational_soc.py",
        "examples/132_recompute_v0250_campaign.py",
        "requirements-pyscf-v250-linux-x86_64-py312.txt",
        "results/v0250_variational_soc_evidence.json",
        "results/v0250_variational_soc_campaign.json",
        "gaussian_dynamics/variational_soc_dynamics_v250.py",
        "gaussian_dynamics/variational_soc_validation_v250.py",
        "gaussian_dynamics/v250_benchmark.py",
    )
    assert all((root / name).is_file() for name in required)

    acceptance = campaign["acceptance"]
    assert acceptance["passed"] is True
    assert acceptance["inherited_gate_count"] == 400
    assert acceptance["validation_gate_count"] == 45
    assert acceptance["core_gate_count"] == 15
    assert acceptance["new_gate_count"] == 60
    assert acceptance["total_gate_count"] == 460
    assert len(acceptance["checks"]) == 460
    assert all(type(value) is bool for value in acceptance["checks"].values())
    assert all(acceptance["checks"].values())

    claims = campaign["claims"]
    assert claims["restricted_single_packet_tdvp_validated"] is True
    assert claims["symmetric_strang_verlet_coupling_validated"] is True
    assert claims["svd_computed_polar_transport_validated"] is True
    assert claims["complete_spinor_soc_propagation_validated"] is True
    assert claims["coordinate_dependent_complex_gauge_covariance_validated"] is True
    assert claims["full_multi_gaussian_tdvp_validated"] is False
    assert claims["adaptive_gaussian_width_tdvp_validated"] is False
    assert claims["plain_verlet_for_general_tdvp_validated"] is False
    assert claims["coordinate_dependent_mass_verlet_validated"] is False
    assert claims["real_pyscf_soc_trajectory_admitted"] is False
    assert claims["general_ab_initio_soc_dynamics_accuracy_validated"] is False

    assert evidence["schema"] == (
        "gnd-symmetric-variational-soc-validation-v0.25.0"
    )
    assert evidence["trajectory_schema"] == (
        "gnd-symmetric-variational-soc-trajectory-v0.25.0"
    )
    assert evidence["audit"]["passed"] is True
    assert len(evidence["audit"]["checks"]) == 45
    assert len(evidence["convergence_receipts"]) == 4
    assert evidence["decisions"]["polar_algorithm"].startswith("SVD")
    assert "implicit midpoint" in evidence["decisions"]["general_tdvp_integrator"]
