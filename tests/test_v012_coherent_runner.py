import numpy as np

from gaussian_dynamics.dynamic_graph_aims import DynamicGraphTBF
from gaussian_dynamics.dynamic_gauge_graph import AnalyticCI2DFrameProvider
from gaussian_dynamics.coherent_lvc_dynamics_v12 import (
    run_coherent_lvc_gaussians,
)
from gaussian_dynamics.electronic_observables import (
    reduced_electronic_density_analytic_ci_diabatic,
)


def _seed():
    return DynamicGraphTBF(
        uid=0,
        state=1,
        q=np.array([0.55,0.45]),
        p=np.array([0.6,0.8]),
        A=1.2*np.eye(2),
        node=("seed",0),
    )


def test_v012_short_run_is_finite_and_spawns():
    provider=AnalyticCI2DFrameProvider(nuclear_mass_au=20.0)

    out=run_coherent_lvc_gaussians(
        [_seed()],
        [1.0+0j],
        provider=provider,
        dt=2e-4,
        steps=25,
        integrator="cayley",
        spawn_action_threshold=1e-7,
        overlap_block=0.99999,
        child_overlap_block=0.99999,
        max_basis=4,
        max_generation=3,
        children_per_event=2,
        minimum_spawn_separation_steps=2,
        position_shifts=(0.0,0.04,-0.04),
        width_scales=(0.7,1.0,1.5),
        store_every=5,
    )

    assert any(e["kind"]=="optimized_spawn" for e in out["events"])
    assert len(out["final_basis"]) >= 2
    assert np.all(np.isfinite(out["final_coefficients"]))

    norms=np.array([r["norm"] for r in out["records"]])
    assert np.max(np.abs(norms-1.0)) < 2e-4

    rho=reduced_electronic_density_analytic_ci_diabatic(
        out["final_coefficients"],
        out["final_basis"],
        normalize=True,
    )
    assert np.isclose(np.trace(rho),1.0)
    assert np.min(np.linalg.eigvalsh(rho)) > -1e-10


def test_v012_cayley_fixed_one_tbf_has_small_norm_drift():
    provider=AnalyticCI2DFrameProvider(nuclear_mass_au=20.0)

    out=run_coherent_lvc_gaussians(
        [_seed()],
        [1.0+0j],
        provider=provider,
        dt=5e-4,
        steps=20,
        integrator="cayley",
        spawn_action_threshold=1e9,
        max_basis=1,
        store_every=1,
    )

    drift=max(abs(r["norm"]-1.0) for r in out["records"])
    assert drift < 2e-6
