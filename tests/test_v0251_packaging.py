import json
from pathlib import Path

from gaussian_dynamics import __version__


def test_v0251_metadata_evidence_and_claim_boundary_are_consistent():
    root = Path(__file__).resolve().parents[1]
    pyproject = (root / "pyproject.toml").read_text(encoding="utf-8")
    citation = (root / "CITATION.cff").read_text(encoding="utf-8")
    readme = (root / "README.md").read_text(encoding="utf-8")
    campaign = json.loads(
        (root / "results/v0251_multigaussian_tdvp_campaign.json").read_text(
            encoding="utf-8"
        )
    )
    evidence = json.loads(
        (root / "results/v0251_multigaussian_tdvp_evidence.json").read_text(
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
        "docs/releases/v0.25.1/V251_RELEASE_NOTES.md",
        "docs/releases/v0.25.1/V251_MULTIGAUSSIAN_TDVP.md",
        "docs/releases/v0.25.1/V251_METRIC_AND_SOLVER.md",
        "docs/releases/v0.25.1/V251_PROGRAM_ARCHITECTURE.md",
        "docs/releases/v0.25.1/V251_ALGORITHM_COMPLEXITY.md",
        "docs/releases/v0.25.1/V251_VALIDATION.md",
        "docs/releases/v0.25.1/V251_BUILD_VALIDATION.md",
        "docs/21_MULTIGAUSSIAN_TDVP.md",
        "examples/133_recompute_v0251_multigaussian_tdvp.py",
        "examples/134_recompute_v0251_campaign.py",
        "requirements-pyscf-v251-linux-x86_64-py312.txt",
        "results/v0251_multigaussian_tdvp_evidence.json",
        "results/v0251_multigaussian_tdvp_campaign.json",
        "gaussian_dynamics/multigaussian_tdvp_v251.py",
        "gaussian_dynamics/multigaussian_tdvp_validation_v251.py",
        "gaussian_dynamics/v251_benchmark.py",
    )
    assert all((root / name).is_file() for name in required)

    acceptance = campaign["acceptance"]
    assert acceptance["passed"] is True
    assert acceptance["inherited_gate_count"] == 460
    assert acceptance["validation_gate_count"] == 55
    assert acceptance["core_gate_count"] == 20
    assert acceptance["new_gate_count"] == 75
    assert acceptance["total_gate_count"] == 535
    assert len(acceptance["checks"]) == 535
    assert all(type(value) is bool for value in acceptance["checks"].values())
    assert all(acceptance["checks"].values())

    claims = campaign["claims"]
    assert claims["frozen_width_multigaussian_tdvp_metric_validated"] is True
    assert claims["implicit_midpoint_nonlinear_solve_validated"] is True
    assert claims["svd_metric_rank_and_compatible_nullspace_validated"] is True
    assert claims["complete_spinor_quadratic_soc_validated"] is True
    assert claims["gaussian_permutation_covariance_validated"] is True
    assert claims["constant_electronic_gauge_covariance_validated"] is True
    assert claims["adaptive_gaussian_width_tdvp_validated"] is False
    assert claims["dynamic_spawning_validated"] is False
    assert claims["dynamic_pruning_validated"] is False
    assert claims["coordinate_dependent_electronic_gauge_covariance_validated"] is False
    assert claims["multidimensional_multigaussian_tdvp_validated"] is False
    assert claims["real_pyscf_soc_trajectory_admitted"] is False
    assert claims["general_ab_initio_soc_dynamics_accuracy_validated"] is False

    assert evidence["schema"] == (
        "gnd-frozen-width-multigaussian-tdvp-validation-v0.25.1"
    )
    assert evidence["trajectory_schema"] == (
        "gnd-frozen-width-multigaussian-tdvp-trajectory-v0.25.1"
    )
    assert evidence["audit"]["passed"] is True
    assert len(evidence["audit"]["checks"]) == 55
    assert len(evidence["convergence_receipts"]) == 4
    assert "fully implicit midpoint" in evidence["decisions"]["integrator"]
    assert "full SVD" in evidence["decisions"]["metric_solver"]
