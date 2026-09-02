from gaussian_dynamics import run_v0213_release_benchmark


def test_v0213_release_benchmark_passes_without_claiming_soc():
    output = run_v0213_release_benchmark()
    assert output["acceptance"]["passed"]
    assert len(output["acceptance"]["checks"]) == 20
    assert output["inherited_v0212_acceptance"]
    assert output["soc"]["physical_hamiltonian_introduced"] is False
    assert output["soc"]["spin_free_mode_permanent"] is True
    assert output["pyscf"]["runtime_validated"] is False
