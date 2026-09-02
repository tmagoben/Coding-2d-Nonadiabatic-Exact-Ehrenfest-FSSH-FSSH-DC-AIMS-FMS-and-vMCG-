from dataclasses import dataclass

from .locality_graph_v16 import (
    conservative_position_overlap_bound,
)


@dataclass(frozen=True)
class LocalSparseCostEstimate:
    local_degree: int
    active_pairs_before: int
    horizon_steps: int
    additional_pair_factorizations: int
    pair_relative_overhead: float
    estimated_nnz_growth: int
    sparse_relative_overhead: float
    electronic_cost_units: float
    electronic_cache_hit: bool
    condition_multiplier: float
    normalized_incremental_cost: float


@dataclass(frozen=True)
class LocalSparseUtilityScore:
    candidate_index: int
    capture_fraction: float
    utility: float
    normalized_incremental_cost: float
    local_degree: int
    estimated_nnz_growth: int
    electronic_cost_units: float
    electronic_cache_hit: bool
    expanded_condition_number: float
    parent_uid: int
    target_state: int
    label: str


def predicted_local_degree(
    candidate_tbf,
    basis,
    *,
    overlap_threshold,
):
    """Cheap degree prediction using the conservative position-overlap bound."""
    degree=0
    for b in basis:
        if (
            conservative_position_overlap_bound(
                candidate_tbf,b
            )
            >=float(overlap_threshold)
        ):
            degree+=1
    return int(degree)


def estimate_local_sparse_incremental_cost(
    candidate_tbf,
    basis,
    *,
    active_offdiagonal_edges,
    overlap_threshold,
    nstate=2,
    horizon_steps=10,
    current_condition=1.0,
    expanded_condition=1.0,
    condition_penalty_weight=0.15,
    electronic_cost_model=None,
    electronic_cost_weight=1.0,
):
    r"""Estimate incremental work for a local/sparse one-TBF addition.

    Pair work uses the candidate's predicted local degree rather than N.

    For each future step the new TBF can create approximately

        degree + 1

    canonical endpoint pairs (including its diagonal) and the same number at the
    midpoint.

    Sparse matrix work is represented by the added block nonzeros.  This is a
    transparent work proxy, not a theorem for sparse-direct factorization fill-in.
    """
    N=len(basis)
    if N<1:
        raise ValueError("basis cannot be empty.")

    degree=predicted_local_degree(
        candidate_tbf,basis,
        overlap_threshold=overlap_threshold,
    )
    h=max(int(horizon_steps),1)
    s=int(nstate)

    active_pairs_before=(
        N+int(active_offdiagonal_edges)
    )
    additional_pairs=2*h*(degree+1)
    current_pairs=2*h*max(active_pairs_before,1)
    pair_relative=additional_pairs/current_pairs

    # One diagonal sxs block plus both orientations of every local offdiagonal block.
    delta_nnz=s*s*(1+2*degree)
    current_nnz=s*s*(
        N+2*int(active_offdiagonal_edges)
    )
    sparse_relative=delta_nnz/max(current_nnz,1)

    cond0=max(float(current_condition),1.0)
    cond1=max(float(expanded_condition),1.0)
    import math
    log_growth=max(
        math.log10(cond1/cond0),
        0.0,
    )
    condition_multiplier=(
        1.0
        +float(condition_penalty_weight)*log_growth
    )

    if electronic_cost_model is None:
        electronic_cost=0.0
        cache_hit=True
    else:
        estimate=electronic_cost_model.estimate(
            candidate_tbf.q
        )
        electronic_cost=float(estimate.cost_units)
        cache_hit=bool(estimate.cache_hit)

    normalized=(
        pair_relative
        +0.5*sparse_relative
        +float(electronic_cost_weight)*electronic_cost
    )*condition_multiplier

    return LocalSparseCostEstimate(
        local_degree=degree,
        active_pairs_before=int(active_pairs_before),
        horizon_steps=h,
        additional_pair_factorizations=
            int(additional_pairs),
        pair_relative_overhead=
            float(pair_relative),
        estimated_nnz_growth=int(delta_nnz),
        sparse_relative_overhead=
            float(sparse_relative),
        electronic_cost_units=
            float(electronic_cost),
        electronic_cache_hit=cache_hit,
        condition_multiplier=
            float(condition_multiplier),
        normalized_incremental_cost=
            float(normalized),
    )


def rank_local_sparse_candidates(
    residual_ranked,
    dynamic_candidates,
    basis,
    *,
    active_offdiagonal_edges,
    overlap_threshold,
    current_condition,
    horizon_steps=10,
    minimum_capture_fraction=0.0,
    minimum_utility=0.0,
    condition_penalty_weight=0.15,
    electronic_cost_model=None,
    electronic_cost_weight=1.0,
):
    """Benefit / local-sparse-cost reranking of a residual-qualified shortlist."""
    out=[]

    for score in residual_ranked:
        capture=float(score.capture_fraction)
        if capture<float(minimum_capture_fraction):
            continue

        item=dynamic_candidates[
            int(score.candidate_index)
        ]
        child=item.candidate.to_tbf(
            uid=-1,
            node_prefix="v16_cost_probe",
        )

        cost=estimate_local_sparse_incremental_cost(
            child,
            basis,
            active_offdiagonal_edges=
                active_offdiagonal_edges,
            overlap_threshold=overlap_threshold,
            horizon_steps=horizon_steps,
            current_condition=current_condition,
            expanded_condition=
                score.expanded_condition_number,
            condition_penalty_weight=
                condition_penalty_weight,
            electronic_cost_model=
                electronic_cost_model,
            electronic_cost_weight=
                electronic_cost_weight,
        )

        utility=capture/max(
            cost.normalized_incremental_cost,
            1e-30,
        )
        if utility<float(minimum_utility):
            continue

        out.append(
            LocalSparseUtilityScore(
                candidate_index=
                    int(score.candidate_index),
                capture_fraction=capture,
                utility=float(utility),
                normalized_incremental_cost=
                    float(cost.normalized_incremental_cost),
                local_degree=int(cost.local_degree),
                estimated_nnz_growth=
                    int(cost.estimated_nnz_growth),
                electronic_cost_units=
                    float(cost.electronic_cost_units),
                electronic_cache_hit=
                    bool(cost.electronic_cache_hit),
                expanded_condition_number=
                    float(score.expanded_condition_number),
                parent_uid=int(score.parent_uid),
                target_state=int(score.target_state),
                label=str(score.label),
            )
        )

    out.sort(
        key=lambda x:(
            -x.utility,
            -x.capture_fraction,
            x.local_degree,
            x.expanded_condition_number,
            x.candidate_index,
        )
    )
    return out
