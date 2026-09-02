from gaussian_dynamics import (
    run_v020_release_benchmark,
)


def test_v020_release_benchmark_passes():
    out=run_v020_release_benchmark()
    assert out["acceptance"]["passed"]
    assert out["canonical"]["average_sparsity_fraction"]>0.75
    assert out["scaling_fit"]["exact_pair_check_exponent"]<1.1
    assert out["controller_demo"]["final_audit_passed"]
