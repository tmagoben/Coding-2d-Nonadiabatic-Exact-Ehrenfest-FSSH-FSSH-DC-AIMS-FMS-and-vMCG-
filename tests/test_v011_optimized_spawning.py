import numpy as np

from gaussian_dynamics.dynamic_graph_aims import DynamicGraphTBF
from gaussian_dynamics.dynamic_gauge_graph import AnalyticCI2DFrameProvider
from gaussian_dynamics.optimized_spawning import (
    classical_energy,
    generate_spawn_candidates,
    select_spawn_children,
)


def parent():
    return DynamicGraphTBF(
        uid=0,
        state=1,
        q=np.array([0.55,0.45]),
        p=np.array([0.6,0.8]),
        A=1.2*np.eye(2),
        node=("seed",0),
    )


def test_spawn_candidates_conserve_classical_energy():
    provider=AnalyticCI2DFrameProvider(nuclear_mass_au=20.0)
    p=parent()

    candidates=generate_spawn_candidates(
        p,
        target=0,
        provider=provider,
        basis=[p],
        position_shifts=(0.0,0.05,-0.05),
        width_scales=(0.7,1.0,1.4),
        overlap_block=0.99999,
    )

    assert candidates

    e0=classical_energy(p.q,p.p,p.state,provider)
    for c in candidates:
        ec=classical_energy(c.q,c.p,c.target_state,provider)
        assert abs(ec-e0) < 1e-10
        assert abs(c.energy_residual) < 1e-10


def test_spawn_search_is_deterministic_and_ranked():
    provider=AnalyticCI2DFrameProvider(nuclear_mass_au=20.0)
    p=parent()

    a=generate_spawn_candidates(
        p,0,provider,[p],
        position_shifts=(0.0,0.05,-0.05),
        width_scales=(0.65,1.0,1.55),
        overlap_block=0.99999,
    )
    b=generate_spawn_candidates(
        p,0,provider,[p],
        position_shifts=(0.0,0.05,-0.05),
        width_scales=(0.65,1.0,1.55),
        overlap_block=0.99999,
    )

    assert len(a)==len(b)
    assert np.allclose([x.score for x in a],[x.score for x in b])
    assert all(a[i].score>=a[i+1].score for i in range(len(a)-1))


def test_multi_child_selection_can_return_distinct_width_or_phase_space_children():
    provider=AnalyticCI2DFrameProvider(nuclear_mass_au=20.0)
    p=parent()

    selected=select_spawn_children(
        p,0,provider,[p],
        children_per_event=2,
        child_overlap_block=0.99999,
        position_shifts=(0.0,0.08,-0.08),
        width_scales=(0.55,1.0,1.8),
        overlap_block=0.999999,
    )

    assert len(selected)>=1
    assert all(c.score>0 for c in selected)

    if len(selected)==2:
        same_q=np.allclose(selected[0].q,selected[1].q)
        same_p=np.allclose(selected[0].p,selected[1].p)
        same_A=np.allclose(selected[0].A,selected[1].A)
        assert not (same_q and same_p and same_A)
