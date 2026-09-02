from gaussian_dynamics import run_v0214_release_benchmark


def test_v0214_release_benchmark_passes_without_claiming_physical_soc():
    output = run_v0214_release_benchmark()

    assert output["acceptance"]["passed"]
    assert len(output["acceptance"]["checks"]) == 21
    assert all(output["acceptance"]["checks"].values())
    assert output["inherited_v0213_acceptance"]
    assert output["soc"]["physical_hamiltonian_introduced"] is False
    assert output["soc"]["physical_derivative_introduced"] is False
    assert output["soc"]["spin_free_mode_permanent"] is True
    assert output["soc"]["first_physical_soc_target"] == "v0.22"
    assert output["pyscf"]["runtime_validated"] is False

