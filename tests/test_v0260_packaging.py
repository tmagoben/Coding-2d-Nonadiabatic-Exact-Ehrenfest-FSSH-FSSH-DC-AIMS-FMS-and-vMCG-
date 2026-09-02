import json
from pathlib import Path

from gaussian_dynamics import __version__


def test_v0260_metadata_artifacts_and_claim_boundaries_are_consistent():
    root = Path(__file__).resolve().parents[1]
    pyproject = (root / "pyproject.toml").read_text(encoding="utf-8")
    citation = (root / "CITATION.cff").read_text(encoding="utf-8")
    readme = (root / "README.md").read_text(encoding="utf-8")
    campaign = json.loads((root / "results/v0260_multidimensional_campaign.json").read_text())
    evidence = json.loads((root / "results/v0260_multidimensional_evidence.json").read_text())
    assert __version__ == "0.27.0"
    assert 'version = "0.27.0"' in pyproject
    assert 'pyscf = ["pyscf==2.13.1"]' in pyproject
    assert "version: 0.27.0" in citation
    assert "date-released: 2026-08-25" in citation
    assert "Current release: v0.27.0" in readme
    required = (
        "V260_RELEASE_NOTES.md",
        "V260_MULTIDIMENSIONAL_TDVP.md",
        "V260_EXACT_CI_SOC_REFERENCE.md",
        "V260_LIFECYCLE_POLICY.md",
        "V260_PROGRAM_ARCHITECTURE.md",
        "V260_ALGORITHM_COMPLEXITY.md",
        "V260_VALIDATION.md",
        "V260_BUILD_VALIDATION.md",
        "docs/24_REFERENCE_FIRST_MULTIDIMENSIONAL_CI_SOC.md",
        "examples/139_recompute_v0260_multidimensional.py",
        "examples/140_recompute_v0260_campaign.py",
        "examples/141_run_v0260_ci_soc_demo.py",
        "requirements-pyscf-v260-linux-x86_64-py312.txt",
        "results/v0260_multidimensional_evidence.json",
        "results/v0260_multidimensional_campaign.json",
        "gaussian_dynamics/multidimensional_soc_v260.py",
        "gaussian_dynamics/multidimensional_gaussian_tdvp_v260.py",
        "gaussian_dynamics/multidimensional_basis_adaptation_v260.py",
        "gaussian_dynamics/multidimensional_validation_v260.py",
        "gaussian_dynamics/v260_benchmark.py",
    )
    assert all((root / name).is_file() for name in required)
    acceptance = campaign["acceptance"]
    assert acceptance["passed"] is True
    assert acceptance["inherited_gate_count"] == 715
    assert acceptance["validation_gate_count"] == 80
    assert acceptance["core_gate_count"] == 30
    assert acceptance["new_gate_count"] == 110
    assert acceptance["total_gate_count"] == 825
    assert len(acceptance["checks"]) == 825
    assert all(type(value) is bool for value in acceptance["checks"].values())
    assert all(acceptance["checks"].values())
    claims = campaign["claims"]
    assert claims["exact_grid"]["two_dimensional_ci_soc_exact_grid_validated"] is True
    assert claims["tdvp"]["multidimensional_adaptive_width_tdvp_validated"] is True
    assert claims["basis"]["multidimensional_residual_driven_spawning_validated"] is True
    assert claims["tdvp"]["full_correlated_width_matrices_validated"] is False
    assert claims["basis"]["full_aims_branching_validated"] is False
    assert claims["basis"]["real_pyscf_soc_trajectory_admitted"] is False
    assert evidence["schema"] == "gnd-multidimensional-ci-soc-validation-v0.26.0"
    assert evidence["passed"] is True
    assert evidence["check_count"] == 80
