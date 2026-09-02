from dataclasses import dataclass
import numpy as np

from .gaussian_nd import gaussian_nd
from .residual_basis_v13 import GaussianCandidate, nuclear_overlap_matrix
from .optimized_spawning import generate_spawn_candidates


@dataclass(frozen=True)
class DynamicDefectCandidate:
    candidate: GaussianCandidate
    parent_uid: int
    target_state: int
    energy_residual: float
    source: str


@dataclass(frozen=True)
class DynamicDefectScore:
    candidate_index: int
    captured_defect_norm: float
    capture_fraction: float
    orthogonal_norm: float
    expanded_condition_number: float
    parent_uid: int
    target_state: int
    label: str


def _candidate_key(c):
    return (
        tuple(np.round(np.asarray(c.q,float),12)),
        tuple(np.round(np.asarray(c.p,float),12)),
        tuple(np.round(np.asarray(c.A,float).reshape(-1),12)),
        int(c.state),
    )


def generate_energy_conserving_defect_candidates(
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
    r"""Generate physically guided candidates, but rank them later by TDSE defect.

    Candidate geometry/momentum construction reuses the v0.11 energy-conserving local
    spawning machinery.  Its coupling-based score is deliberately ignored here.

    For every parent TBF the dictionary may include:

    - same-guidance-surface candidates: basis/dispersion enrichment;
    - other-guidance-surface candidates: branching-like trajectory enrichment.

    The spinor-complete quantum basis attached to every candidate still contains both
    electronic amplitudes.  `target_state` controls only the future classical guidance
    force of that nuclear Gaussian center.
    """
    if len(basis)==0:
        return []

    nstate=len(provider.evaluate(np.asarray(basis[0].q,float)).energies)
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
                        f"pos={item.position_direction}:{item.position_shift:g};"
                        f"mom={item.momentum_direction};"
                        f"width={item.width_scale:g}"
                    ),
                )
                key=_candidate_key(candidate)
                if key in seen:
                    continue
                seen.add(key)

                out.append(
                    DynamicDefectCandidate(
                        candidate=candidate,
                        parent_uid=int(parent.uid),
                        target_state=int(target),
                        energy_residual=float(item.energy_residual),
                        source="energy_conserving_local_dictionary",
                    )
                )

    return out


def rank_dynamic_defect_candidates_prepared(
    defect,
    basis,
    dynamic_candidates,
    grid,
    condition_limit=1e8,
    orthogonal_norm_floor=1e-8,
    exact_condition_top=32,
    max_return=8,
):
    r"""Vectorized TDSE-defect ranking for a dynamic candidate dictionary.

    Let B be the N x G matrix of current nuclear Gaussian values and Q be the K x G
    candidate matrix.  The implementation evaluates all candidate residual overlaps
    and all candidate projections against the current span with dense matrix products:

        X = B^dag Q
        alpha = S^-1 X
        n_perp = ||q||^2 - x^dag alpha
        b_perp = <q|R> - alpha^dag <B|R>.

    The residual-capture score is

        Delta_k = ||b_perp,k||^2 / n_perp,k.

    Exact expanded overlap condition numbers are evaluated only for the strongest
    residual candidates, not for all K entries.
    """
    if len(dynamic_candidates)==0:
        return []

    area=float(grid.area)
    points=grid.points
    G=int(np.prod(points.shape[:-1]))

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

    S=(B.conj()@B.T)*area
    S=0.5*(S+S.conj().T)
    X=(B.conj()@Q.T)*area

    alpha=np.linalg.lstsq(S,X,rcond=1e-12)[0]

    qnorm=np.real(np.sum(np.abs(Q)**2,axis=1)*area)
    nperp=qnorm-np.real(
        np.sum(np.conj(X)*alpha,axis=0)
    )
    nperp=np.maximum(nperp,0.0)

    R=np.asarray(defect.residual,dtype=complex).reshape(G,-1)
    Bres=(B.conj()@R)*area
    Qres=(Q.conj()@R)*area
    bperp=Qres-alpha.conj().T@Bres

    captured=np.full(len(dynamic_candidates),-np.inf)
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

    # Exact analytic nuclear overlap conditioning is checked for the strongest
    # candidates only.  If none of that shortlist is admissible, the remaining
    # candidates are checked in descending score order as a safe fallback.
    shortlist=list(order[:max(1,int(exact_condition_top))])
    shortlist += [
        int(i) for i in order[max(1,int(exact_condition_top)):]
    ]

    Snuc=nuclear_overlap_matrix(basis)
    scores=[]

    for idx in shortlist:
        item=dynamic_candidates[int(idx)]
        c=item.candidate

        s=np.array([
            # Analytic Gaussian overlap for final conditioning decision.
            # Candidate electronic guidance state is irrelevant to nuclear S.
            _analytic_overlap(b,c)
            for b in basis
        ],dtype=complex)

        Sexp=np.empty(
            (len(basis)+1,len(basis)+1),
            dtype=complex,
        )
        Sexp[:-1,:-1]=Snuc
        Sexp[:-1,-1]=s
        Sexp[-1,:-1]=np.conj(s)
        Sexp[-1,-1]=1.0
        cond=float(np.linalg.cond(Sexp))
        if not np.isfinite(cond) or cond>float(condition_limit):
            continue

        value=float(max(captured[idx],0.0))
        scores.append(
            DynamicDefectScore(
                candidate_index=int(idx),
                captured_defect_norm=float(np.sqrt(value)),
                capture_fraction=float(
                    value/max(defect.residual_norm**2,1e-30)
                ),
                orthogonal_norm=float(nperp[idx]),
                expanded_condition_number=cond,
                parent_uid=int(item.parent_uid),
                target_state=int(item.target_state),
                label=str(c.label),
            )
        )
        if len(scores)>=max(1,int(max_return)):
            break

    scores.sort(
        key=lambda x:(-x.capture_fraction,x.candidate_index)
    )
    return scores


def _analytic_overlap(basis_tbf, candidate):
    from .gaussian_general import gaussian_overlap_general

    return gaussian_overlap_general(
        basis_tbf.q,basis_tbf.p,basis_tbf.A,
        candidate.q,candidate.p,candidate.A,
    )
