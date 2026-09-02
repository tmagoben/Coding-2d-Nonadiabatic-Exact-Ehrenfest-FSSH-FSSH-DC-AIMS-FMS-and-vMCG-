import numpy as np

from gaussian_dynamics.analytic_molecular_backend_v19 import (
    AnalyticMolecularLVCBackendV19,
    AnalyticMolecularLVCConfigV19,
    default_diatomic_two_mode_map_v19,
)
from gaussian_dynamics.benchmark_provider_nd import (
    LVC2DGeneralizedProvider,
)
from gaussian_dynamics.molecular_direct_provider_v19 import (
    TrackedMolecularDirectProviderV19,
    BackendEvaluationPolicyV19,
)


def test_molecular_cartesian_projection_matches_analytic_lvc():
    gmap=default_diatomic_two_mode_map_v19()
    backend=AnalyticMolecularLVCBackendV19(gmap)
    provider=TrackedMolecularDirectProviderV19(
        backend,gmap
    )
    mass=float(
        backend.generalized_mass_matrix_au[0,0]
    )
    ref=LVC2DGeneralizedProvider(
        nuclear_mass_au=mass
    )

    for q in (
        np.array([-0.7,0.3]),
        np.array([-0.2,0.4]),
        np.array([0.4,0.35]),
    ):
        a=provider.evaluate(q)
        b=ref.evaluate(q)
        assert np.allclose(a.energies,b.energies,atol=1e-12)
        assert np.allclose(a.gradients_q,b.gradients_q,atol=1e-12)
        assert np.allclose(a.nac_q,b.nac_q,atol=1e-12)
        assert np.allclose(
            a.mass_matrix_q_au,
            b.mass_matrix_q_au,
            atol=1e-10,
        )


def test_tracking_removes_deliberate_root_swaps_and_sign_flips():
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

    path=[
        np.array([-0.8,0.35]),
        np.array([-0.4,0.35]),
        np.array([-0.1,0.35]),
        np.array([0.2,0.35]),
        np.array([0.5,0.35]),
    ]

    for q in path:
        a=clean.evaluate(q)
        b=scrambled.evaluate(q)
        assert np.allclose(a.energies,b.energies,atol=1e-12)
        assert np.allclose(a.gradients_q,b.gradients_q,atol=1e-12)
        assert np.allclose(a.nac_q,b.nac_q,atol=1e-12)

    assert scrambled.diagnostics.tracking_ambiguities==0


def test_cache_and_cost_estimate_are_explicit():
    gmap=default_diatomic_two_mode_map_v19()
    backend=AnalyticMolecularLVCBackendV19(gmap)
    provider=TrackedMolecularDirectProviderV19(
        backend,gmap
    )

    q=np.array([-0.4,0.3])
    provider.evaluate(q)
    calls=backend.calls
    provider.evaluate(q.copy())

    assert backend.calls==calls
    assert provider.diagnostics.cache_hits==1
    assert provider.cost_estimate(q)["cache_hit"]
    assert (
        provider.cost_estimate(
            q+np.array([0.01,0.0])
        )["nearby_cache"]
    )


def test_bounded_nearest_cache_failure_fallback_is_opt_in():
    gmap=default_diatomic_two_mode_map_v19()
    backend=AnalyticMolecularLVCBackendV19(
        gmap,
        AnalyticMolecularLVCConfigV19(
            fail_if_q0_greater_than=0.1
        ),
    )
    provider=TrackedMolecularDirectProviderV19(
        backend,gmap,
        failure=BackendEvaluationPolicyV19(
            failure_policy="nearest_cache",
            max_fallback_distance=0.05,
        ),
    )

    provider.evaluate(np.array([0.08,0.3]))
    p=provider.evaluate(np.array([0.12,0.3]))

    assert p.metadata["v19_failure_fallback"]
    assert p.metadata["fallback_distance"]<0.05
    assert provider.diagnostics.fallback_uses==1


def test_nearest_anchor_tracking_is_order_tolerant_after_reference_seed():
    gmap=default_diatomic_two_mode_map_v19()
    backend=AnalyticMolecularLVCBackendV19(
        gmap,
        AnalyticMolecularLVCConfigV19(
            scramble_roots=True
        ),
    )
    provider=TrackedMolecularDirectProviderV19(
        backend,gmap
    )
    mass=float(
        backend.generalized_mass_matrix_au[0,0]
    )
    ref=LVC2DGeneralizedProvider(
        nuclear_mass_au=mass
    )

    points=[
        np.array([x,0.35])
        for x in np.linspace(-0.8,0.8,17)
    ]
    provider.evaluate(points[0])
    order=[
        8,4,12,2,14,6,10,1,
        15,3,13,5,11,7,9,16,
    ]

    for i in order:
        a=provider.evaluate(points[i])
        b=ref.evaluate(points[i])
        assert np.allclose(
            a.energies,b.energies,atol=1e-12
        )
        assert np.allclose(
            a.nac_q,b.nac_q,atol=1e-11
        )

    assert provider.diagnostics.tracking_ambiguities==0
