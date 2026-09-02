import json
from pathlib import Path

from gaussian_dynamics import __version__


def test_v0252_metadata_evidence_and_claim_boundary_are_consistent():
    root = Path(__file__).resolve().parents[1]
    pyproject = (root / "pyproject.toml").read_text(encoding="utf-8")
    citation = (root / "CITATION.cff").read_text(encoding="utf-8")
    readme = (root / "README.md").read_text(encoding="utf-8")
    campaign = json.loads(
        (root / "results/v0252_adaptive_multigaussian_tdvp_campaign.json").read_text(
            encoding="utf-8"
        )
    )
    evidence = json.loads(
        (root / "results/v0252_adaptive_multigaussian_tdvp_evidence.json").read_text(
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
        "V252_RELEASE_NOTES.md",
        "V252_ADAPTIVE_MULTIGAUSSIAN_TDVP.md",
        "V252_WIDTH_AND_SOLVER_POLICY.md",
        "V252_PROGRAM_ARCHITECTURE.md",
        "V252_ALGORITHM_COMPLEXITY.md",
        "V252_VALIDATION.md",
        "V252_BUILD_VALIDATION.md",
        "docs/22_ADAPTIVE_MULTIGAUSSIAN_TDVP.md",
        "examples/135_recompute_v0252_adaptive_tdvp.py",
        "examples/136_recompute_v0252_campaign.py",
        "requirements-pyscf-v252-linux-x86_64-py312.txt",
        "results/v0252_adaptive_multigaussian_tdvp_evidence.json",
        "results/v0252_adaptive_multigaussian_tdvp_campaign.json",
        "gaussian_dynamics/adaptive_multigaussian_tdvp_v252.py",
        "gaussian_dynamics/adaptive_multigaussian_tdvp_validation_v252.py",
        "gaussian_dynamics/v252_benchmark.py",
    )
    assert all((root / name).is_file() for name in required)

    acceptance = campaign["acceptance"]
    assert acceptance["passed"] is True
    assert acceptance["inherited_gate_count"] == 535
    assert acceptance["validation_gate_count"] == 70
    assert acceptance["core_gate_count"] == 25
    assert acceptance["new_gate_count"] == 95
    assert acceptance["total_gate_count"] == 630
    assert len(acceptance["checks"]) == 630
    assert all(type(value) is bool for value in acceptance["checks"].values())
    assert all(acceptance["checks"].values())

    claims = campaign["claims"]
    assert claims["adaptive_width_multigaussian_tdvp_validated"] is True
    assert claims["log_width_positivity_and_quadratic_chirp_validated"] is True
    assert claims["implicit_midpoint_adaptive_nonlinear_solve_validated"] is True
    assert claims["single_packet_thawed_harmonic_reduction_validated"] is True
    assert claims["frozen_coherent_state_reduction_validated"] is True
    assert claims["dynamic_spawning_validated"] is False
    assert claims["dynamic_pruning_validated"] is False
    assert claims["coordinate_dependent_electronic_gauge_covariance_validated"] is False
    assert claims["multidimensional_adaptive_width_tdvp_validated"] is False
    assert claims["full_correlated_width_matrices_validated"] is False
    assert claims["real_pyscf_soc_trajectory_admitted"] is False
    assert claims["general_ab_initio_soc_dynamics_accuracy_validated"] is False

    assert evidence["schema"] == (
        "gnd-adaptive-width-multigaussian-tdvp-validation-v0.25.2"
    )
    assert evidence["trajectory_schema"] == (
        "gnd-adaptive-width-multigaussian-tdvp-trajectory-v0.25.2"
    )
    assert evidence["audit"]["passed"] is True
    assert len(evidence["audit"]["checks"]) == 70
    assert len(evidence["convergence_receipts"]) == 4
    assert evidence["decisions"]["analytic_moment_degree"] == 4
    assert "fully implicit midpoint" in evidence["decisions"]["integrator"]
    assert "full SVD" in evidence["decisions"]["metric_solver"]
