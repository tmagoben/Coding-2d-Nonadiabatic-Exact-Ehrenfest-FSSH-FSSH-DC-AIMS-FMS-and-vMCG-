from dataclasses import dataclass
import numpy as np

from .gaussian_nd import gaussian_nd
from .residual_basis_v13 import GaussianCandidate
from .optimized_spawning import generate_spawn_candidates
from .pair_cache_v15 import GaussianPairCache


@dataclass(frozen=True)
class DynamicDefectCandidateV15:
    candidate: GaussianCandidate
    parent_uid: int
    target_state: int
    energy_residual: float
    source: str


@dataclass
class CachedDynamicDefectScore:
    candidate_index: int
    captured_defect_norm: float
    capture_fraction: float
    orthogonal_norm: float
    expanded_condition_number: float
    parent_uid: int
    target_state: int
    label: str
    expanded_cache: GaussianPairCache


def _candidate_key(c):
    return (
        tuple(np.round(np.asarray(c.q,float),12)),
        tuple(np.round(np.asarray(c.p,float),12)),
        tuple(np.round(np.asarray(c.A,float).reshape(-1),12)),
        int(c.state),
    )


def generate_energy_conserving_defect_candidates_v15(
    basis,
    provider,
    position_shifts=(0.0,0.06,-0.06),
    width_scales=(0.75,1.0,1.35),
    momentum_directions=("nac","momentum"),
    include_same_surface=True,
    include_other_surfaces=True,
    overlap_block=0.999999,
    energy_tolerance=1e-10,
):
    """v0.15 local physical candidate dictionary.

    Placement remains energy-conserving/locally guided.  Selection is deferred to the
    measured TDSE-defect and cost-aware utility.
    """
    if len(basis)==0:
        return []

    nstate=len(
        provider.evaluate(
            np.asarray(basis[0].q,float)
        ).energies
    )
    out=[]
    seen=set()

    for parent in basis:
        targets=[]
        if include_same_surface:
            targets.append(int(parent.state))
        if include_other_surfaces:
            targets.extend(
                s for s in range(nstate)
                if s!=int(parent.state)
            )

        for target in targets:
            raw=generate_spawn_candidates(
                parent,
                target,
                provider,
                basis,
                position_shifts=position_shifts,
                width_scales=width_scales,
                momentum_directions=momentum_directions,
                overlap_block=overlap_block,
                novelty_power=0.0,
                energy_tolerance=energy_tolerance,
            )

            for item in raw:
                candidate=GaussianCandidate(
                    q=np.asarray(item.q,float),
                    p=np.asarray(item.p,float),
                    A=np.asarray(item.A,float),
                    state=int(target),
                    label=(
                        f"parent={int(parent.uid)};"
                        f"target={int(target)};"
                        f"pos={item.position_direction}:"
                        f"{item.position_shift:g};"
                        f"mom={item.momentum_direction};"
                        f"width={item.width_scale:g}"
                    ),
                )
                key=_candidate_key(candidate)
                if key in seen:
                    continue
                seen.add(key)

                out.append(
                    DynamicDefectCandidateV15(
                        candidate=candidate,
                        parent_uid=int(parent.uid),
                        target_state=int(target),
                        energy_residual=float(item.energy_residual),
                        source="energy_conserving_local_dictionary",
                    )
                )

    return out


def rank_dynamic_defect_candidates_cached(
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
):
    r"""Residual-rank candidates while reusing the endpoint Gaussian-pair snapshot.

    The expensive K x G residual contractions are still vectorized.

    Exact overlap conditioning is evaluated only for a residual shortlist.  For each
    shortlisted candidate, a temporary expanded pair cache inherits every old-old pair
    from the current endpoint snapshot and computes only the N+1 new child pairs.
    If that candidate is accepted, the same expanded cache is returned and reused by
    the incremental matrix expansion.
    """
    if len(dynamic_candidates)==0:
        return []

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

    Q=np.asarray([
        gaussian_nd(
            points,
            item.candidate.q,
            item.candidate.p,
            item.candidate.A,
        ).reshape(G)
        for item in dynamic_candidates
    ],dtype=complex)

    Sgrid=(B.conj()@B.T)*area
    Sgrid=0.5*(Sgrid+Sgrid.conj().T)
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

    R=np.asarray(
        defect.residual,dtype=complex
    ).reshape(G,-1)
    Bres=(B.conj()@R)*area
    Qres=(Q.conj()@R)*area
    bperp=Qres-alpha.conj().T@Bres

    captured=np.full(
        len(dynamic_candidates),
        -np.inf,
    )
    valid=nperp>=float(orthogonal_norm_floor)
    captured[valid]=(
        np.sum(np.abs(bperp[valid])**2,axis=1)
        /nperp[valid]
    )

    finite=np.flatnonzero(np.isfinite(captured))
    if len(finite)==0:
        return []

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
            node_prefix="v15_candidate",
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
                orthogonal_norm=float(nperp[idx]),
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
    return scores[:max(1,int(max_return))]
