import numpy as np
from scipy.optimize import linear_sum_assignment

from .state_tracking import StateTrackingResult


def _best_assignment(score,forbidden=None):
    """Maximum-weight perfect matching; `forbidden=(i,j)` excludes one edge."""
    W=np.asarray(score,dtype=float).copy()
    if W.ndim!=2 or W.shape[0]!=W.shape[1]:
        raise ValueError("score matrix must be square.")
    n=W.shape[0]

    cost=-W
    if forbidden is not None:
        i,j=forbidden
        # A finite penalty avoids optimizer issues with infinities.
        scale=max(float(np.max(np.abs(cost))),1.0)
        cost[int(i),int(j)]=scale*(n+1)*1e6

    rows,cols=linear_sum_assignment(cost)
    if not np.array_equal(rows,np.arange(n)):
        order=np.argsort(rows)
        cols=cols[order]
    value=float(np.sum(W[np.arange(n),cols]))

    if forbidden is not None:
        i,j=forbidden
        if int(cols[int(i)])==int(j):
            return None,-np.inf
    return np.asarray(cols,dtype=int),value


def scalable_maximum_overlap_assignment_v19(
    overlap,
    minimum_overlap=0.50,
    minimum_score_margin=0.05,
    real_gauge=True,
    imaginary_tolerance=1e-8,
):
    """Exact maximum-overlap assignment without n! permutation enumeration.

    The best assignment maximizes sum_i |O_i,perm_i|^2 by the Hungarian algorithm.

    The exact second-best score is obtained by forbidding each edge used by the best
    assignment in turn and taking the best remaining perfect matching. Any distinct
    assignment must omit at least one edge of the best assignment, so the maximum over
    these constrained optima is the true second-best score.

    Complexity:
        best matching: O(n^3)
        n constrained matchings: O(n^4)
    """
    O=np.asarray(overlap,dtype=complex)
    if O.ndim!=2 or O.shape[0]!=O.shape[1]:
        raise ValueError("State-overlap matrix must be square.")
    n=O.shape[0]
    W=np.abs(O)**2

    perm,best_score=_best_assignment(W)
    assigned=O[np.arange(n),perm].copy()

    if n<=1:
        second_score=-np.inf
    else:
        alternatives=[]
        for i,j in enumerate(perm):
            _,score=_best_assignment(
                W,forbidden=(i,int(j))
            )
            alternatives.append(score)
        second_score=float(max(alternatives))

    margin=(
        np.inf
        if not np.isfinite(second_score)
        else float(best_score-second_score)
    )

    reasons=[]
    mags=np.abs(assigned)
    if np.min(mags)<float(minimum_overlap):
        reasons.append(
            f"minimum assigned overlap {np.min(mags):.6g} "
            f"is below threshold {minimum_overlap:.6g}"
        )
    if np.isfinite(margin) and margin<float(minimum_score_margin):
        reasons.append(
            f"assignment score margin {margin:.6g} "
            f"is below threshold {minimum_score_margin:.6g}"
        )

    phases=np.ones(n,dtype=complex)
    for i,z in enumerate(assigned):
        if abs(z)<1e-15:
            continue
        if real_gauge:
            if abs(np.imag(z))>imaginary_tolerance:
                reasons.append(
                    f"assigned overlap for tracked state {i} has imaginary part "
                    f"{np.imag(z):.6g}, incompatible with real-gauge tracking"
                )
            phases[i]=1.0 if np.real(z)>=0.0 else -1.0
        else:
            phases[i]=np.conj(z)/abs(z)

    assigned_after=assigned*phases

    return StateTrackingResult(
        permutation=perm,
        phase_factors=phases,
        assigned_overlaps=assigned_after,
        best_score=float(best_score),
        second_best_score=float(second_score),
        score_margin=float(margin),
        ambiguous=bool(reasons),
        reasons=tuple(reasons),
    )
