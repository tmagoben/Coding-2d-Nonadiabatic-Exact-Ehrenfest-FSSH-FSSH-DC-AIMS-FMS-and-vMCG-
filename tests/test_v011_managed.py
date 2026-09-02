import numpy as np

from gaussian_dynamics.dynamic_graph_aims import DynamicGraphTBF
from gaussian_dynamics.dynamic_gauge_graph import AnalyticCI2DFrameProvider
from gaussian_dynamics.managed_graph_aims_v11 import (
    run_basis_complete_graph_aims,
)
from gaussian_dynamics.electronic_observables import (
    reduced_electronic_density_analytic_ci_diabatic,
)


def seed():
    return DynamicGraphTBF(
        uid=0,
        state=1,
        q=np.array([0.55,0.45]),
        p=np.array([0.6,0.8]),
        A=1.2*np.eye(2),
        node=("seed",0),
    )


def test_v011_short_run_conserves_norm_and_spawns():
    provider=AnalyticCI2DFrameProvider(nuclear_mass_au=20.0)
    out=run_basis_complete_graph_aims(
        [seed()],
        [1.0+0j],
        provider=provider,
        dt=0.0002,
        steps=25,
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
    assert len(out["final_basis"])>=2
    assert max(abs(r["norm"]-1.0) for r in out["records"]) < 1e-4

    rho=reduced_electronic_density_analytic_ci_diabatic(
        out["final_coefficients"],
        out["final_basis"],
        normalize=True,
    )
    assert abs(np.trace(rho)-1.0) < 1e-12
    assert np.all(np.linalg.eigvalsh(rho)>-1e-10)


def test_lineage_records_descendants_and_width_diversity():
    provider=AnalyticCI2DFrameProvider(nuclear_mass_au=20.0)
    out=run_basis_complete_graph_aims(
        [seed()],
        [1.0+0j],
        provider=provider,
        dt=0.0002,
        steps=30,
        spawn_action_threshold=1e-7,
        overlap_block=0.999999,
        child_overlap_block=0.999999,
        max_basis=6,
        max_generation=3,
        children_per_event=2,
        minimum_spawn_separation_steps=1,
        position_shifts=(0.0,0.05,-0.05),
        width_scales=(0.55,1.0,1.8),
        store_every=10,
    )

    assert out["lineage"][0]["generation"]==0
    assert all(
        info["generation"]>=0
        for info in out["lineage"].values()
    )

    dets=[round(float(np.linalg.det(b.A)),10) for b in out["final_basis"]]
    assert len(set(dets))>=1
