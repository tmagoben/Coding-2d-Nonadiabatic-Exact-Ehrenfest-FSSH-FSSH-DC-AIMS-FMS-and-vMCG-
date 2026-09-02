from dataclasses import dataclass
import numpy as np

from .gaussian_nd import gaussian_nd
from .defect_candidates_v15 import CachedDynamicDefectScore


@dataclass(frozen=True)
class BatchedRankingDiagnosticsV18:
    candidate_count: int
    grid_points: int
    batch_size: int
    batches: int
    basis_grid_elements: int
    peak_candidate_grid_elements: int
    dense_candidate_grid_elements: int
    peak_grid_element_reduction_fraction: float

    def as_dict(self):
        return {
            "candidate_count":int(self.candidate_count),
            "grid_points":int(self.grid_points),
            "batch_size":int(self.batch_size),
            "batches":int(self.batches),
            "basis_grid_elements":
                int(self.basis_grid_elements),
            "peak_candidate_grid_elements":
                int(self.peak_candidate_grid_elements),
            "dense_candidate_grid_elements":
                int(self.dense_candidate_grid_elements),
            "peak_grid_element_reduction_fraction":
                float(self.peak_grid_element_reduction_fraction),
        }


def rank_dynamic_defect_candidates_batched_v18(
    defect,
    basis,
    dynamic_candidates,
    grid,
    current_cache,
    Snuc,
    *,
    condition_limit=1e8,
    orthogonal_norm_floor=1e-8,
    exact_condition_top=12,
    max_return=8,
    batch_size=32,
    return_diagnostics=False,
):
    r"""Exact v0.15 residual-ranking algebra with bounded candidate-grid memory.

    v0.15 forms all candidate Gaussian fields as a K x G complex array.

    v0.18 retains the same projection equations but processes candidates in batches of
    B, so peak candidate-grid storage becomes O(B G) instead of O(K G).

    Only scalar capture/orthogonality values for all K candidates are retained before
    exact conditioning is applied to the residual shortlist.
    """
    K=len(dynamic_candidates)
    if K==0:
        empty=[]
        if return_diagnostics:
            return empty,BatchedRankingDiagnosticsV18(
                candidate_count=0,
                grid_points=int(np.prod(grid.points.shape[:-1])),
                batch_size=max(int(batch_size),1),
                batches=0,
                basis_grid_elements=0,
                peak_candidate_grid_elements=0,
                dense_candidate_grid_elements=0,
                peak_grid_element_reduction_fraction=0.0,
            )
        return empty

    batch_size=max(int(batch_size),1)
    area=float(grid.area)
    points=grid.points
    G=int(np.prod(points.shape[:-1]))
    n=len(basis)

    B=np.asarray([
        gaussian_nd(
            points,b.q,b.p,b.A
        ).reshape(G)
        for b in basis
    ],dtype=complex)

    Sgrid=(B.conj()@B.T)*area
    Sgrid=0.5*(Sgrid+Sgrid.conj().T)

    R=np.asarray(
        defect.residual,dtype=complex
    ).reshape(G,-1)
    Bres=(B.conj()@R)*area

    captured=np.full(K,-np.inf,dtype=float)
    nperp_all=np.zeros(K,dtype=float)

    batches=0
    peak_candidate_elements=0

    for start in range(0,K,batch_size):
        stop=min(start+batch_size,K)
        items=dynamic_candidates[start:stop]
        Q=np.asarray([
            gaussian_nd(
                points,
                item.candidate.q,
                item.candidate.p,
                item.candidate.A,
            ).reshape(G)
            for item in items
        ],dtype=complex)

        batches+=1
        peak_candidate_elements=max(
            peak_candidate_elements,
            int(Q.size),
        )

        X=(B.conj()@Q.T)*area
        alpha=np.linalg.lstsq(
            Sgrid,X,rcond=1e-12
        )[0]

        qnorm=np.real(
            np.sum(np.abs(Q)**2,axis=1)*area
        )
        nperp=qnorm-np.real(
            np.sum(np.conj(X)*alpha,axis=0)
        )
        nperp=np.maximum(nperp,0.0)

        Qres=(Q.conj()@R)*area
        bperp=Qres-alpha.conj().T@Bres

        valid=nperp>=float(orthogonal_norm_floor)
        local=np.full(len(items),-np.inf,dtype=float)
        local[valid]=(
            np.sum(np.abs(bperp[valid])**2,axis=1)
            /nperp[valid]
        )

        captured[start:stop]=local
        nperp_all[start:stop]=nperp

    finite=np.flatnonzero(np.isfinite(captured))
    if len(finite)==0:
        scores=[]
    else:
        order=finite[
            np.argsort(captured[finite])[::-1]
        ]
        order=order[:max(
            int(exact_condition_top),
            int(max_return),
            1,
        )]

        Snuc=np.asarray(Snuc,dtype=complex)
        scores=[]

        for idx in order:
            item=dynamic_candidates[int(idx)]
            child=item.candidate.to_tbf(
                uid=-1,
                node_prefix="v18_candidate",
            )

            expanded_cache=current_cache.expanded(
                child
            )
            s=np.array([
                expanded_cache.pair(i,n).overlap
                for i in range(n)
            ],dtype=complex)
            diag=expanded_cache.pair(n,n).overlap

            Sexp=np.empty(
                (n+1,n+1),
                dtype=complex,
            )
            Sexp[:-1,:-1]=Snuc
            Sexp[:-1,-1]=s
            Sexp[-1,:-1]=np.conj(s)
            Sexp[-1,-1]=diag

            cond=float(np.linalg.cond(Sexp))
            if (
                not np.isfinite(cond)
                or cond>float(condition_limit)
            ):
                continue

            value=float(max(captured[idx],0.0))
            scores.append(
                CachedDynamicDefectScore(
                    candidate_index=int(idx),
                    captured_defect_norm=float(
                        np.sqrt(value)
                    ),
                    capture_fraction=float(
                        value/max(
                            defect.residual_norm**2,
                            1e-30,
                        )
                    ),
                    orthogonal_norm=float(
                        nperp_all[idx]
                    ),
                    expanded_condition_number=cond,
                    parent_uid=int(item.parent_uid),
                    target_state=int(item.target_state),
                    label=str(item.candidate.label),
                    expanded_cache=expanded_cache,
                )
            )

        scores.sort(
            key=lambda x:(
                -x.capture_fraction,
                x.expanded_condition_number,
                x.candidate_index,
            )
        )
        scores=scores[:max(1,int(max_return))]

    dense_elements=int(K*G)
    reduction=(
        1.0-peak_candidate_elements/max(dense_elements,1)
    )
    diag=BatchedRankingDiagnosticsV18(
        candidate_count=K,
        grid_points=G,
        batch_size=batch_size,
        batches=batches,
        basis_grid_elements=int(B.size),
        peak_candidate_grid_elements=
            int(peak_candidate_elements),
        dense_candidate_grid_elements=
            dense_elements,
        peak_grid_element_reduction_fraction=
            float(reduction),
    )

    if return_diagnostics:
        return scores,diag
    return scores
