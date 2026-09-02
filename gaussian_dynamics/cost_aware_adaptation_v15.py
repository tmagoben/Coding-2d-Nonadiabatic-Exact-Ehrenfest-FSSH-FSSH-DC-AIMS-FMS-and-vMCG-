from dataclasses import dataclass
import math


@dataclass(frozen=True)
class IncrementalCostEstimate:
    n_basis_before: int
    n_basis_after: int
    horizon_steps: int
    defect_checks: int

    additional_pair_factorizations: int
    additional_cayley_cubic_units: int
    additional_defect_cubic_units: int
    pair_relative_overhead: float
    solve_relative_overhead: float
    defect_relative_overhead: float

    condition_multiplier: float
    normalized_incremental_cost: float
    estimated_seconds: float | None


@dataclass(frozen=True)
class CostAwareCandidateScore:
    candidate_index: int
    capture_fraction: float
    expanded_condition_number: float
    normalized_incremental_cost: float
    estimated_incremental_seconds: float | None
    utility: float
    parent_uid: int
    target_state: int
    label: str


def estimate_one_tbf_incremental_cost(
    n_basis,
    *,
    nstate=2,
    horizon_steps=10,
    defect_checks=1,
    current_condition=1.0,
    expanded_condition=1.0,
    condition_penalty_weight=0.15,
    pair_seconds_per_factorization=None,
    cayley_seconds_per_cubic_unit=None,
):
    r"""Estimate cost of carrying one additional TBF for a short control horizon.

    The dense reference runner builds, per propagation step:

      - one NEW endpoint pair snapshot,
      - one midpoint pair snapshot,
      - one Cayley solve.

    Adding one TBF changes the canonical pair count by N+1 for each snapshot.

    It also changes dense solve work from

        (sN)^3

    to

        [s(N+1)]^3.

    A defect checkpoint contributes one additional projected dense solve.  Grid
    reconstruction overhead is represented by the same O(1/N) relative-defect term
    rather than converted into seconds unless a dedicated timing model is available.

    `condition_multiplier` is a numerical-stability penalty, not literal CPU time.
    """
    N=int(n_basis)
    s=int(nstate)
    h=max(int(horizon_steps),1)
    q=max(int(defect_checks),0)

    if N<1:
        raise ValueError("n_basis must be positive.")

    additional_pairs=2*h*(N+1)

    current_pair_work=(
        2*h*(N*(N+1)//2)
    )
    pair_relative=additional_pairs/max(
        current_pair_work,1
    )

    m0=s*N
    m1=s*(N+1)
    delta_cube=m1**3-m0**3

    additional_cayley=h*delta_cube
    current_cayley=h*(m0**3)
    solve_relative=additional_cayley/max(
        current_cayley,1
    )

    additional_defect=q*delta_cube
    current_defect=q*(m0**3)
    defect_relative=(
        additional_defect/max(current_defect,1)
        if q>0 else 0.0
    )

    cond0=max(float(current_condition),1.0)
    cond1=max(float(expanded_condition),1.0)
    log_growth=max(
        math.log10(cond1/cond0),
        0.0,
    )
    condition_multiplier=(
        1.0
        +float(condition_penalty_weight)*log_growth
    )

    normalized=(
        pair_relative
        +solve_relative
        +0.25*defect_relative
    )*condition_multiplier

    estimated_seconds=None
    if (
        pair_seconds_per_factorization is not None
        and cayley_seconds_per_cubic_unit is not None
    ):
        estimated_seconds=(
            additional_pairs
            *max(float(pair_seconds_per_factorization),0.0)
            +(
                additional_cayley
                +additional_defect
            )
            *max(float(cayley_seconds_per_cubic_unit),0.0)
        )*condition_multiplier

    return IncrementalCostEstimate(
        n_basis_before=N,
        n_basis_after=N+1,
        horizon_steps=h,
        defect_checks=q,
        additional_pair_factorizations=
            int(additional_pairs),
        additional_cayley_cubic_units=
            int(additional_cayley),
        additional_defect_cubic_units=
            int(additional_defect),
        pair_relative_overhead=float(pair_relative),
        solve_relative_overhead=float(solve_relative),
        defect_relative_overhead=float(defect_relative),
        condition_multiplier=float(condition_multiplier),
        normalized_incremental_cost=float(normalized),
        estimated_seconds=(
            None if estimated_seconds is None
            else float(estimated_seconds)
        ),
    )


def rank_candidates_by_cost_aware_utility(
    residual_ranked,
    *,
    n_basis,
    current_condition,
    horizon_steps,
    defect_checks=1,
    minimum_capture_fraction=0.0,
    minimum_utility=0.0,
    condition_penalty_weight=0.15,
    pair_seconds_per_factorization=None,
    cayley_seconds_per_cubic_unit=None,
):
    r"""Rerank an already residual-qualified short list by benefit / incremental cost.

    Benefit is the predicted fraction of the current squared TDSE defect captured by
    the candidate.

    Cost includes:

      - additional canonical pair factorizations over the control horizon;
      - increase in dense Cayley/Galerkin cubic work;
      - one or more defect-check dense solves;
      - a logarithmic penalty for worsening overlap conditioning.

    No exact future-time quantum observable is used.
    """
    out=[]

    for score in residual_ranked:
        capture=float(score.capture_fraction)
        if capture<float(minimum_capture_fraction):
            continue

        estimate=estimate_one_tbf_incremental_cost(
            n_basis,
            horizon_steps=horizon_steps,
            defect_checks=defect_checks,
            current_condition=current_condition,
            expanded_condition=
                score.expanded_condition_number,
            condition_penalty_weight=
                condition_penalty_weight,
            pair_seconds_per_factorization=
                pair_seconds_per_factorization,
            cayley_seconds_per_cubic_unit=
                cayley_seconds_per_cubic_unit,
        )

        denominator=max(
            estimate.normalized_incremental_cost,
            1e-30,
        )
        utility=capture/denominator

        if utility<float(minimum_utility):
            continue

        out.append(
            CostAwareCandidateScore(
                candidate_index=int(score.candidate_index),
                capture_fraction=capture,
                expanded_condition_number=float(
                    score.expanded_condition_number
                ),
                normalized_incremental_cost=float(
                    estimate.normalized_incremental_cost
                ),
                estimated_incremental_seconds=
                    estimate.estimated_seconds,
                utility=float(utility),
                parent_uid=int(score.parent_uid),
                target_state=int(score.target_state),
                label=str(score.label),
            )
        )

    out.sort(
        key=lambda x:(
            -x.utility,
            -x.capture_fraction,
            x.expanded_condition_number,
            x.candidate_index,
        )
    )
    return out
