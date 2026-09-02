from gaussian_dynamics import run_v0220_release_benchmark


def test_v0220_release_benchmark_passes_with_physical_analytic_soc_only():
    output = run_v0220_release_benchmark()

    assert output["acceptance"]["passed"]
    assert output["acceptance"]["inherited_gate_count"] == 21
    assert output["acceptance"]["new_gate_count"] == 32
    assert output["acceptance"]["total_gate_count"] == 53
    assert len(output["acceptance"]["checks"]) == 53
    assert all(output["acceptance"]["checks"].values())
    assert output["inherited_v0214"]["acceptance"]["passed"]
    assert output["soc"]["physical_hamiltonian_introduced"] is True
    assert output["soc"]["physical_derivative_introduced"] is True
    assert output["soc"]["analytic_models_only"] is True
    assert output["soc"]["ab_initio_SOC_validated"] is False
    assert output["pyscf"]["runtime_validated"] is False

