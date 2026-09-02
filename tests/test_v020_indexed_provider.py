import numpy as np

from gaussian_dynamics.indexed_molecular_provider_v20 import (
    BufferedKDTreeIndexV20,
    IndexedTrackedMolecularDirectProviderV20,
)
from gaussian_dynamics.analytic_molecular_backend_v19 import (
    AnalyticMolecularLVCBackendV19,
    default_diatomic_two_mode_map_v19,
)


def test_buffered_kdtree_nearest_matches_bruteforce():
    rng=np.random.default_rng(20020)
    idx=BufferedKDTreeIndexV20(
        dimension=3,rebuild_batch=5
    )
    points={}
    for i in range(23):
        q=rng.normal(size=3)
        points[i]=q
        idx.insert(i,q)

    for _ in range(50):
        q=rng.normal(size=3)
        key,d=idx.nearest(q)
        brute=min(
            points,
            key=lambda k:(
                np.linalg.norm(q-points[k]),
                repr(k),
            ),
        )
        db=float(np.linalg.norm(q-points[brute]))
        assert key==brute
        assert abs(d-db)<1e-12

    assert idx.diagnostics.kd_queries>0
    assert idx.diagnostics.rebuilds>=1


def test_indexed_molecular_provider_matches_v19_tracking():
    gmap=default_diatomic_two_mode_map_v19()
    provider=IndexedTrackedMolecularDirectProviderV20(
        AnalyticMolecularLVCBackendV19(gmap),
        gmap,
        rebuild_batch=4,
    )

    points=[
        np.array([x,0.35])
        for x in np.linspace(-0.8,0.8,17)
    ]
    provider.evaluate(points[0])
    for i in [8,4,12,2,14,6,10,1,15,3,13,5,11,7,9,16]:
        p=provider.evaluate(points[i])
        assert np.all(np.isfinite(p.energies))

    diag=provider.diagnostics_dict()["spatial_index"]
    assert diag["nearest_queries"]>=16
    assert diag["kd_queries"]>=1
