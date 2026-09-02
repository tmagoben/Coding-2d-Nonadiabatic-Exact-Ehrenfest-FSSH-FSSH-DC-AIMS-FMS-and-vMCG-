from gaussian_dynamics.convergence_campaign_v18 import (
    observed_order_from_dt,
    axis_sensitivity_summary,
    refinement_ladder_summary,
    successive_self_convergence_order,
)


def test_observed_order_recovers_second_order_example():
    rows=[
        {"coordinates":{"dt":0.02},"err":4e-4},
        {"coordinates":{"dt":0.01},"err":1e-4},
        {"coordinates":{"dt":0.005},"err":2.5e-5},
    ]
    out=observed_order_from_dt(rows,"err")
    assert abs(out[0]["observed_order"]-2.0)<1e-12
    assert abs(out[1]["observed_order"]-2.0)<1e-12


def test_axis_and_refinement_summaries():
    rows=[
        {"err":0.3},
        {"err":0.2},
        {"err":0.1},
    ]
    axis=axis_sensitivity_summary(rows,"err",higher_is_better=False)
    assert abs(axis["span"]-0.2)<1e-15
    assert axis["best_index"]==2

    fidelity=[
        {"f":0.8},
        {"f":0.9},
        {"f":0.85},
    ]
    fsum=axis_sensitivity_summary(
        fidelity,"f",higher_is_better=True
    )
    assert fsum["best_index"]==1
    assert fsum["worst_index"]==0

    ref=refinement_ladder_summary(rows,"err")
    assert ref["fine_better_than_coarse"]
    assert ref["strictly_monotone"]


def test_successive_self_convergence_order():
    # For second order and r=2, successive differences differ by a factor of four.
    p=successive_self_convergence_order(4e-4,1e-4,2.0)
    assert abs(p-2.0)<1e-12
