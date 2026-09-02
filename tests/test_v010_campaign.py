import numpy as np

from gaussian_dynamics.benchmark_campaign import (
    CIPassageConfig,
    run_exact_grid_timestep_surface,
    select_finest_exact_reference,
    run_managed_parameter_surface,
)


def tiny_config():
    return CIPassageConfig(
        q0=(0.55,0.45),
        p0=(0.6,0.8),
        A_diag=(1.2,1.2),
        mass=20.0,
        final_time=0.002,
        half_width=3.0,
    )


def test_exact_surface_and_reference_selection():
    rows=run_exact_grid_timestep_surface(
        tiny_config(),
        grid_values=(16,20),
        dt_values=(0.001,0.0005),
    )
    assert len(rows)==4
    ref=select_finest_exact_reference(rows)
    assert ref["grid_n"]==20
    assert ref["dt"]==0.0005
    assert abs(np.sum(ref["populations"])-1.0) < 1e-10


def test_managed_parameter_surface_has_expected_cartesian_product():
    rows=run_managed_parameter_surface(
        tiny_config(),
        dts=(0.001,),
        spa_orders=(0,1),
        spawn_action_thresholds=(1e-6,),
        max_basis_values=(2,3),
    )
    assert len(rows)==4
    assert {r["spa_order"] for r in rows}=={0,1}
    assert {r["max_basis"] for r in rows}=={2,3}


def test_repeated_spawning_can_grow_basis_beyond_one_child():
    config=CIPassageConfig(final_time=0.18)
    from gaussian_dynamics.benchmark_campaign import run_managed_passage

    out=run_managed_passage(
        config,
        dt=0.005,
        max_basis=4,
        overlap_block=0.9999,
        minimum_spawn_separation_steps=5,
        store_every=18,
    )

    spawn_count=sum(e["kind"]=="spawn" for e in out["events"])
    assert spawn_count >= 2
    assert len(out["final_basis"]) >= 3
