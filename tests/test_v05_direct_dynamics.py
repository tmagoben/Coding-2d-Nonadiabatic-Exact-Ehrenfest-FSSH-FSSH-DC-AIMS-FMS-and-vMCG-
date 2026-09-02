import numpy as np

from gaussian_dynamics.local_gaussian_nd import LocalAdiabaticTBF
from gaussian_dynamics.benchmark_provider_nd import LVC2DGeneralizedProvider
from gaussian_dynamics.direct_dynamics_nd import (
    energy_conserving_child,
    maybe_spawn_once,
    run_backend_spawned_gaussians,
)


def test_general_mass_energy_conserving_child():
    provider=LVC2DGeneralizedProvider(nuclear_mass_au=20.0)
    parent=LocalAdiabaticTBF(
        state=1,
        q=np.array([0.55,0.45]),
        p=np.array([0.6,0.8]),
        A=np.eye(2),
    )

    child=energy_conserving_child(parent,0,provider)
    assert child is not None

    point=provider.evaluate(parent.q)
    B=np.linalg.inv(point.mass_matrix_q_au)

    before=0.5*parent.p@B@parent.p+point.energies[parent.state]
    after=0.5*child.p@B@child.p+point.energies[child.state]

    assert abs(before-after) < 1e-12


def test_short_backend_driven_spawned_run():
    provider=LVC2DGeneralizedProvider(nuclear_mass_au=20.0)
    parent=LocalAdiabaticTBF(
        state=1,
        q=np.array([0.55,0.45]),
        p=np.array([0.6,0.8]),
        A=1.2*np.eye(2),
    )

    out=run_backend_spawned_gaussians(
        [parent],
        [1.0+0j],
        provider,
        dt=0.0002,
        steps=30,
        spawn_threshold=1e-6,
        overlap_block=0.9,
        max_basis=2,
        store_every=5,
    )

    assert out["basis_size"][-1] == 2
    assert len(out["events"]) == 1
    assert np.all(np.isfinite(out["norm"]))
    assert np.max(np.abs(out["norm"]-1.0)) < 2e-5
    assert out["state_populations"][-1,0] > 0.0
