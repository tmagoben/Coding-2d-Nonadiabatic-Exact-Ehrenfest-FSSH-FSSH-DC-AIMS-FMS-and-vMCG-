from dataclasses import dataclass

from gaussian_dynamics.cost_aware_adaptation_v15 import (
    estimate_one_tbf_incremental_cost,
    rank_candidates_by_cost_aware_utility,
)


@dataclass(frozen=True)
class _Score:
    candidate_index: int
    capture_fraction: float
    expanded_condition_number: float
    parent_uid: int
    target_state: int
    label: str


def test_incremental_cost_contains_pair_and_cubic_growth():
    c=estimate_one_tbf_incremental_cost(
        10,
        horizon_steps=10,
        defect_checks=1,
        current_condition=100.0,
        expanded_condition=200.0,
    )

    assert c.additional_pair_factorizations==220
    assert c.additional_cayley_cubic_units>0
    assert c.additional_defect_cubic_units>0
    assert c.normalized_incremental_cost>0.0
    assert c.condition_multiplier>1.0


def test_cost_aware_ranking_penalizes_bad_conditioning():
    scores=[
        _Score(
            candidate_index=0,
            capture_fraction=0.15,
            expanded_condition_number=1e5,
            parent_uid=0,target_state=1,
            label="high-capture-bad-condition",
        ),
        _Score(
            candidate_index=1,
            capture_fraction=0.145,
            expanded_condition_number=150.0,
            parent_uid=0,target_state=1,
            label="slightly-lower-capture-stable",
        ),
    ]

    ranked=rank_candidates_by_cost_aware_utility(
        scores,
        n_basis=10,
        current_condition=100.0,
        horizon_steps=10,
        defect_checks=1,
        condition_penalty_weight=1.0,
    )

    assert ranked[0].candidate_index==1
    assert ranked[0].utility>ranked[1].utility


def test_cost_gate_can_reject_low_benefit_candidate():
    scores=[
        _Score(
            candidate_index=0,
            capture_fraction=1e-5,
            expanded_condition_number=100.0,
            parent_uid=0,target_state=1,
            label="tiny-benefit",
        )
    ]

    ranked=rank_candidates_by_cost_aware_utility(
        scores,
        n_basis=10,
        current_condition=100.0,
        horizon_steps=10,
        defect_checks=1,
        minimum_utility=0.01,
    )
    assert ranked==[]
