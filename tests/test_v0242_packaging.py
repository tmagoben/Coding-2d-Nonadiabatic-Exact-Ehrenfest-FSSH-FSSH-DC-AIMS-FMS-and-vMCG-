import json
from pathlib import Path

from gaussian_dynamics import __version__


def test_v0242_release_metadata_evidence_and_claim_boundary_are_consistent():
    root = Path(__file__).resolve().parents[1]
    pyproject = (root / "pyproject.toml").read_text(encoding="utf-8")
    citation = (root / "CITATION.cff").read_text(encoding="utf-8")
    readme = (root / "README.md").read_text(encoding="utf-8")
    campaign = json.loads(
        (root / "results/v0242_pyscf_differential_soc_campaign.json").read_text(
            encoding="utf-8"
        )
    )
    evidence = json.loads(
        (root / "results/v0242_pyscf_differential_soc_evidence.json").read_text(
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
        "docs/releases/v0.24.2/V242_RELEASE_NOTES.md",
        "docs/releases/v0.24.2/V242_PYSCF_DIFFERENTIAL_SOC.md",
        "docs/releases/v0.24.2/V242_PROGRAM_ARCHITECTURE.md",
        "docs/releases/v0.24.2/V242_ALGORITHM_COMPLEXITY.md",
        "docs/releases/v0.24.2/V242_VALIDATION.md",
        "docs/releases/v0.24.2/V242_BUILD_VALIDATION.md",
        "docs/19_PYSCF_CONNECTED_GEOMETRY_SOC.md",
        "examples/129_recompute_v0242_pyscf_differential_soc.py",
        "examples/130_recompute_v0242_campaign.py",
        "requirements-pyscf-v242-linux-x86_64-py312.txt",
        "results/v0242_pyscf_differential_soc_evidence.json",
        "results/v0242_pyscf_differential_soc_campaign.json",
        "gaussian_dynamics/pyscf_differential_soc_v242.py",
        "gaussian_dynamics/v242_benchmark.py",
    )
    assert all((root / name).is_file() for name in required)

    acceptance = campaign["acceptance"]
    assert acceptance["passed"] is True
    assert acceptance["inherited_gate_count"] == 315
    assert acceptance["runtime_gate_count"] == 60
    assert acceptance["core_gate_count"] == 25
    assert acceptance["new_gate_count"] == 85
    assert acceptance["total_gate_count"] == 400
    assert len(acceptance["checks"]) == 400
    assert all(type(value) is bool for value in acceptance["checks"].values())
    assert all(acceptance["checks"].values())

    claims = campaign["claims"]
    assert claims["direct_jk_somf_execution_validated"] is True
    assert claims["rank_five_tensor_avoided_in_production_path"] is True
    assert claims["connected_geometry_soc_snapshots_validated"] is True
    assert claims["complete_doublet_overlap_transport_validated"] is True
    assert claims["transported_spin_free_derivative_preview_validated"] is True
    assert claims["transported_soc_derivative_preview_validated"] is True
    assert claims["continuous_physical_derivative_connection_validated"] is False
    assert claims["full_cartesian_derivative_tensor_validated"] is False
    assert claims["analytic_soc_derivatives_validated"] is False
    assert claims["real_mixed_multiplicity_runtime_validated"] is False
    assert claims["trajectory_ready_molecular_soc_validated"] is False
    assert claims["live_molecular_soc_backend_admitted"] is False
    assert claims["ab_initio_soc_accuracy_validated"] is False

    assert evidence["schema"] == (
        "gnd-pyscf-connected-geometry-soc-differential-v0.24.2"
    )
    assert evidence["audit"]["passed"] is True
    assert len(evidence["audit"]["checks"]) == 60
    assert evidence["audit"]["metrics"]["retained_endpoint_snapshot_count"] == 6
    assert len(evidence["scan"]["endpoint_snapshots"]) == 3
    assert len(evidence["scan"]["derivative_records"]) == 3
    assert evidence["scan"]["capability"] == (
        "connected_geometry_differential_preview"
    )
    assert evidence["scan"]["direct_jk_explicit_max_abs_error"] < 2.0e-14
    assert evidence["scan"]["convergence_metrics"]["K_soc"][
        "finest_norm_frobenius"
    ] > 1.0e-4
