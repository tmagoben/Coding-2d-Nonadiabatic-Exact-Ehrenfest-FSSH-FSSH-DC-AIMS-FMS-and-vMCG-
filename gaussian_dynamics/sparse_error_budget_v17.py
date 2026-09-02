import numpy as np

from .edge_importance_v17 import (
    EdgeImportanceSettingsV17,
    ErrorControlledGaussianLocalityGraphV17,
)
from .sparse_pair_matrices_v16 import (
    build_sparse_spinor_lvc_matrices,
    audit_sparse_lvc_matrices_against_dense,
)


def monotone_nonincreasing(values,atol=1e-12):
    x=np.asarray(values,dtype=float)
    if x.size<2:
        return True
    return bool(np.all(x[1:]<=x[:-1]+float(atol)))


def score_threshold_snapshot_sweep(
    basis,
    provider,
    *,
    dt,
    enter_scores=(0.12,0.08,0.06,0.04,0.03,0.02,0.01),
    search_overlap_floor=1e-5,
    overlap_weight=1.0,
    hamiltonian_weight=0.20,
    time_connection_weight=1.0,
):
    """Dense-audited snapshot convergence as the S/H/T edge threshold is relaxed.

    The local L2 importance budget is disabled here so the sweep isolates the edge
    score itself.
    """
    rows=[]

    for enter in enter_scores:
        enter=float(enter)
        graph=ErrorControlledGaussianLocalityGraphV17(
            provider,
            dt,
            EdgeImportanceSettingsV17(
                enter_score=enter,
                exit_score=0.5*enter,
                search_overlap_floor=
                    search_overlap_floor,
                overlap_weight=overlap_weight,
                hamiltonian_weight=
                    hamiltonian_weight,
                time_connection_weight=
                    time_connection_weight,
                local_omitted_score_l2_budget=
                    1e30,
            ),
        )
        update=graph.update(basis)
        mats=build_sparse_spinor_lvc_matrices(
            update,provider
        )
        audit=audit_sparse_lvc_matrices_against_dense(
            basis,provider,mats
        )

        rows.append({
            "enter_score":enter,
            "active_edges":
                int(update.active_offdiagonal_edges),
            "edge_fraction":
                float(update.edge_fraction),
            "omitted_score_l2":
                float(update.omitted_candidate_score_l2),
            "omitted_score_max":
                float(update.omitted_candidate_score_max),
            **audit,
        })

    return rows


def local_score_budget_snapshot_sweep(
    basis,
    provider,
    *,
    dt,
    enter_score=0.06,
    budgets=(1e30,0.10,0.08,0.05,0.03,0.01),
    search_overlap_floor=1e-5,
    overlap_weight=1.0,
    hamiltonian_weight=0.20,
    time_connection_weight=1.0,
):
    """Dense-audited snapshot sweep of the global local-importance L2 budget proxy."""
    rows=[]

    for budget in budgets:
        graph=ErrorControlledGaussianLocalityGraphV17(
            provider,
            dt,
            EdgeImportanceSettingsV17(
                enter_score=float(enter_score),
                exit_score=0.5*float(enter_score),
                search_overlap_floor=
                    search_overlap_floor,
                overlap_weight=overlap_weight,
                hamiltonian_weight=
                    hamiltonian_weight,
                time_connection_weight=
                    time_connection_weight,
                local_omitted_score_l2_budget=
                    float(budget),
            ),
        )
        update=graph.update(basis)
        mats=build_sparse_spinor_lvc_matrices(
            update,provider
        )
        audit=audit_sparse_lvc_matrices_against_dense(
            basis,provider,mats
        )

        rows.append({
            "budget":float(budget),
            "active_edges":
                int(update.active_offdiagonal_edges),
            "budget_promoted_edges":
                int(update.budget_promoted_edges),
            "omitted_score_l2":
                float(update.omitted_candidate_score_l2),
            **audit,
        })

    return rows


def summarize_snapshot_convergence(
    threshold_rows,
    budget_rows,
):
    threshold_H=[
        row["relative_H_frobenius_error"]
        for row in threshold_rows
    ]
    threshold_S=[
        row["relative_S_frobenius_error"]
        for row in threshold_rows
    ]

    # Budget rows are ordered from loose -> strict. Error should fall or stay equal.
    budget_H=[
        row["relative_H_frobenius_error"]
        for row in budget_rows
    ]
    budget_S=[
        row["relative_S_frobenius_error"]
        for row in budget_rows
    ]

    return {
        "threshold_S_monotone":
            monotone_nonincreasing(threshold_S),
        "threshold_H_monotone":
            monotone_nonincreasing(threshold_H),
        "budget_S_monotone":
            monotone_nonincreasing(budget_S),
        "budget_H_monotone":
            monotone_nonincreasing(budget_H),
        "finest_threshold_S_error":
            float(threshold_S[-1]),
        "finest_threshold_H_error":
            float(threshold_H[-1]),
        "strictest_budget_S_error":
            float(budget_S[-1]),
        "strictest_budget_H_error":
            float(budget_H[-1]),
    }
