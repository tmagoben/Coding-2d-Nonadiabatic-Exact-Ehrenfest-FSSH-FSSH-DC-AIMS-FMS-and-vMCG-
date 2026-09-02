import json
from pathlib import Path

from gaussian_dynamics import __version__


def test_v0253_metadata_evidence_and_claim_boundary_are_consistent():
    root = Path(__file__).resolve().parents[1]
    pyproject = (root / "pyproject.toml").read_text(encoding="utf-8")
    citation = (root / "CITATION.cff").read_text(encoding="utf-8")
    readme = (root / "README.md").read_text(encoding="utf-8")
    campaign = json.loads(
        (root / "results/v0253_controlled_basis_campaign.json").read_text(encoding="utf-8")
    )
    evidence = json.loads(
        (root / "results/v0253_controlled_basis_evidence.json").read_text(encoding="utf-8")
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
        "V253_RELEASE_NOTES.md",
        "V253_CONTROLLED_BASIS_ADAPTATION.md",
        "V253_LIFECYCLE_POLICY.md",
        "V253_PROGRAM_ARCHITECTURE.md",
        "V253_ALGORITHM_COMPLEXITY.md",
        "V253_VALIDATION.md",
        "V253_BUILD_VALIDATION.md",
        "docs/23_CONTROLLED_BASIS_ADAPTATION.md",
        "examples/137_recompute_v0253_controlled_basis.py",
        "examples/138_recompute_v0253_campaign.py",
        "requirements-pyscf-v253-linux-x86_64-py312.txt",
        "results/v0253_controlled_basis_evidence.json",
        "results/v0253_controlled_basis_campaign.json",
        "gaussian_dynamics/controlled_basis_adaptation_v253.py",
        "gaussian_dynamics/controlled_basis_validation_v253.py",
        "gaussian_dynamics/v253_benchmark.py",
    )
    assert all((root / name).is_file() for name in required)
    acceptance = campaign["acceptance"]
    assert acceptance["passed"] is True
    assert acceptance["inherited_gate_count"] == 630
    assert acceptance["validation_gate_count"] == 60
    assert acceptance["core_gate_count"] == 25
    assert acceptance["new_gate_count"] == 85
    assert acceptance["total_gate_count"] == 715
    assert len(acceptance["checks"]) == 715
    assert all(type(value) is bool for value in acceptance["checks"].values())
    assert all(acceptance["checks"].values())
    claims = campaign["claims"]
    assert claims["controlled_residual_driven_spawning_validated"] is True
    assert claims["coefficient_only_newborn_activation_validated"] is True
    assert claims["projection_guarded_pruning_validated"] is True
    assert claims["overlap_projection_guarded_merging_validated"] is True
    assert claims["general_aims_branching_validated"] is False
    assert claims["multidimensional_spawning_validated"] is False
    assert claims["coordinate_dependent_electronic_gauge_covariance_validated"] is False
    assert claims["real_pyscf_soc_trajectory_admitted"] is False
    assert claims["general_ab_initio_soc_dynamics_accuracy_validated"] is False
    assert evidence["schema"] == "gnd-controlled-basis-validation-v0.25.3"
    assert evidence["trajectory_schema"] == "gnd-controlled-basis-trajectory-v0.25.3"
    assert evidence["audit"]["passed"] is True
    assert len(evidence["audit"]["checks"]) == 60
    assert "full-SVD" in evidence["decisions"]["projection_policy"]
