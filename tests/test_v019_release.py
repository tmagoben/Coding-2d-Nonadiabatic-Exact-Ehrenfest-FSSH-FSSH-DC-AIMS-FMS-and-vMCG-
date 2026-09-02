from gaussian_dynamics import run_v019_release_benchmark


def test_v019_release_benchmark_passes():
    out=run_v019_release_benchmark()
    assert out["acceptance"]["passed"]
    assert out["acceptance"]["checks"]["order_tolerant_tracking"]
    assert out["state_tracking_scaling"]["nstate"]==16
