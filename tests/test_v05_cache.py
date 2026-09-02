import numpy as np

from gaussian_dynamics.backend_cache import DiskCachedGeneralizedProvider
from gaussian_dynamics.benchmark_provider_nd import LVC2DGeneralizedProvider


def test_disk_cache_reuses_point(tmp_path):
    base=LVC2DGeneralizedProvider()
    cached=DiskCachedGeneralizedProvider(base,tmp_path,namespace="test")

    q=np.array([0.7,0.4])
    a=cached.evaluate(q)
    b=cached.evaluate(q)

    assert cached.misses == 1
    assert cached.hits == 1
    assert np.allclose(a.energies,b.energies)
    assert np.allclose(a.gradients_q,b.gradients_q)
    assert np.allclose(a.nac_q,b.nac_q)
