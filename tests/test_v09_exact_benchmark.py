import numpy as np

from gaussian_dynamics.exact_benchmark import run_exact_ci_reference


def test_exact_reference_adiabatic_populations_sum_to_norm():
    out=run_exact_ci_reference(grid_n=28,dt=0.002,final_time=0.004)
    pops=out['final_populations_adiabatic']
    assert abs(np.sum(pops)-out['norm']) < 2e-10
    assert abs(out['norm']-1.0) < 2e-10
