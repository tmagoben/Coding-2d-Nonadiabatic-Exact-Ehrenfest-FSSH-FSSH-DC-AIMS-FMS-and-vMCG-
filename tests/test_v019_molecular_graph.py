import numpy as np

from gaussian_dynamics.analytic_molecular_backend_v19 import (
    AnalyticMolecularLVCBackendV19,
    AnalyticMolecularLVCConfigV19,
    default_diatomic_two_mode_map_v19,
)
from gaussian_dynamics.molecular_direct_provider_v19 import (
    TrackedMolecularDirectProviderV19,
)
from gaussian_dynamics.local_gaussian_nd import LocalAdiabaticTBF
from gaussian_dynamics.molecular_gauge_graph_v19 import (
    build_molecular_centroid_graph_v19,
)


def _basis():
    A=1.1*np.eye(2)
    return [
        LocalAdiabaticTBF(
            0,np.array([-0.6,0.35]),
            np.array([0.3,0.0]),A
        ),
        LocalAdiabaticTBF(
            1,np.array([0.1,0.40]),
            np.array([0.1,0.1]),A
        ),
        LocalAdiabaticTBF(
            0,np.array([0.7,0.32]),
            np.array([-0.2,0.0]),A
        ),
    ]


def test_center_centroid_graph_is_gauge_stable_under_raw_scrambling():
    gmap=default_diatomic_two_mode_map_v19()

    clean=TrackedMolecularDirectProviderV19(
        AnalyticMolecularLVCBackendV19(gmap),
        gmap,
    )
    scrambled=TrackedMolecularDirectProviderV19(
        AnalyticMolecularLVCBackendV19(
            gmap,
            AnalyticMolecularLVCConfigV19(
                scramble_roots=True
            ),
        ),
        gmap,
    )

    g1=build_molecular_centroid_graph_v19(
        _basis(),clean
    )
    g2=build_molecular_centroid_graph_v19(
        _basis(),scrambled
    )

    S1,H1=g1.matrices()
    S2,H2=g2.matrices()

    assert np.allclose(S1,S2,atol=1e-11)
    assert np.allclose(H1,H2,atol=1e-11)
    assert np.linalg.norm(S1-S1.conj().T)<1e-10
    assert np.linalg.norm(H1-H1.conj().T)<1e-10
