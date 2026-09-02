import numpy as np

from gaussian_dynamics.ensemble_benchmark import ensemble_statistics


def test_ensemble_statistics_mean_std_sem():
    x=np.array([
        [0.2,0.8],
        [0.4,0.6],
        [0.3,0.7],
    ])
    s=ensemble_statistics(x)

    assert np.allclose(s.mean,[0.3,0.7])
    assert s.nsamples==3
    assert np.allclose(s.sem,s.std/np.sqrt(3))
