from gaussian_dynamics import run_v0212_release_benchmark


def test_v0212_release_benchmark_passes():
    out=run_v0212_release_benchmark()
    assert out["acceptance"]["passed"]
    assert out["inherited_v021_acceptance"]
    assert out["pyscf"]["runtime_validated"] is False
    assert out["theme"].startswith("pre-SOC")
