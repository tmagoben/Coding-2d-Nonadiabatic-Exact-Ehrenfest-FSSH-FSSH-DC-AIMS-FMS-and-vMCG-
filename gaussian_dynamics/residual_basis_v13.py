from dataclasses import dataclass
import itertools
import numpy as np

from .dynamic_graph_aims import DynamicGraphTBF
from .gaussian_nd import gaussian_nd
from .gaussian_general import gaussian_overlap_general
from .initial_projection_v12 import (
    project_grid_wavefunction_to_spinor_complete_basis,
)
from .electronic_observables import (
    exact_reduced_electronic_density_diabatic,
)


@dataclass(frozen=True)
class GaussianCandidate:
    q: np.ndarray
    p: np.ndarray
    A: np.ndarray
    state: int
    label: str = ""

    def to_tbf(self, uid, node_prefix="residual"):
        return DynamicGraphTBF(
            uid=int(uid),
            state=int(self.state),
            q=np.asarray(self.q,float).copy(),
            p=np.asarray(self.p,float).copy(),
            A=np.asarray(self.A,float).copy(),
            node=(node_prefix,int(uid),self.label),
        )


@dataclass(frozen=True)
class ResidualCandidateScore:
    candidate_index: int
    predicted_gain: float
    orthogonal_norm: float


@dataclass(frozen=True)
class ResidualSelectionStep:
    basis_size: int
    selected_candidate_index: int
    selected_label: str
    predicted_gain: float
    actual_residual_reduction: float
    projection_fidelity: float
    relative_residual: float
    density_error: float
    condition_number: float


@dataclass
class ResidualBasisBuild:
    basis: list
    projection: object
    history: list
    target_density: np.ndarray
    candidate_indices: list

    @property
    def final_relative_residual(self):
        return float(self.projection.relative_residual)

    @property
    def final_fidelity(self):
        return float(self.projection.fidelity)


def cartesian_offsets_2d(radius=1.0, spacing=0.2):
    """Deterministic square lattice of 2D center offsets."""
    radius=float(radius)
    spacing=float(spacing)
    if radius<0.0 or spacing<=0.0:
        raise ValueError("radius must be nonnegative and spacing positive.")

    n=int(round(radius/spacing))
    vals=np.arange(-n,n+1,dtype=float)*spacing
    return [
        np.array([x,y],dtype=float)
        for x,y in itertools.product(vals,vals)
    ]


def generate_gaussian_dictionary(
    q0,
    p0,
    A0,
    state,
    position_offsets,
    width_scales=(1.0,1.5,2.0,3.0,4.0,6.0),
    momentum_offsets=None,
):
    """Generate a deterministic Gaussian candidate dictionary.

    The dictionary is deliberately explicit rather than optimized by a hidden
    nonlinear solver.  It can therefore be reproduced exactly from the documented
    center offsets, momentum offsets, and width scales.
    """
    q0=np.asarray(q0,float)
    p0=np.asarray(p0,float)
    A0=np.asarray(A0,float)

    if q0.ndim!=1 or p0.shape!=q0.shape:
        raise ValueError("q0 and p0 must be equal-length vectors.")
    if A0.shape!=(len(q0),len(q0)):
        raise ValueError("A0 has incompatible shape.")

    if momentum_offsets is None:
        momentum_offsets=[np.zeros_like(p0)]

    candidates=[]
    seen=set()

    for scale in width_scales:
        scale=float(scale)
        if scale<=0.0:
            raise ValueError("width scales must be positive.")
        A=scale*A0

        for dq in position_offsets:
            dq=np.asarray(dq,float)
            if dq.shape!=q0.shape:
                raise ValueError("position offset has incompatible shape.")

            for dp in momentum_offsets:
                dp=np.asarray(dp,float)
                if dp.shape!=p0.shape:
                    raise ValueError("momentum offset has incompatible shape.")

                q=q0+dq
                p=p0+dp
                key=(
                    tuple(np.round(q,12)),
                    tuple(np.round(p,12)),
                    tuple(np.round(A.reshape(-1),12)),
                    int(state),
                )
                if key in seen:
                    continue
                seen.add(key)

                dq_label=tuple(float(x) for x in np.round(dq,6))
                dp_label=tuple(float(x) for x in np.round(dp,6))
                label=(
                    f"dq={dq_label};"
                    f"dp={dp_label};"
                    f"width_scale={scale:g}"
                )
                candidates.append(
                    GaussianCandidate(
                        q=q.copy(),
                        p=p.copy(),
                        A=A.copy(),
                        state=int(state),
                        label=label,
                    )
                )

    return candidates


def nuclear_overlap_matrix(basis):
    n=len(basis)
    S=np.zeros((n,n),dtype=complex)

    for i,bi in enumerate(basis):
        for j,bj in enumerate(basis):
            S[i,j]=gaussian_overlap_general(
                bi.q,bi.p,bi.A,
                bj.q,bj.p,bj.A,
            )
    return 0.5*(S+S.conj().T)


def normalized_grid_density(psi,dx):
    rho=exact_reduced_electronic_density_diabatic(
        np.asarray(psi,dtype=complex),
        float(dx),
        float(dx),
    )
    tr=np.trace(rho)
    if abs(tr)<1e-15:
        raise ValueError("wavefunction has zero reduced-density trace.")
    return rho/tr


def candidate_orthogonal_norm(
    candidate,
    basis,
    overlap_matrix=None,
    eigenvalue_floor=1e-12,
):
    r"""Norm of the candidate after projection out of the current nuclear span.

    If
        s_i = <g_i|g_c>,
    then
        ||g_c^\perp||^2 = 1 - s^\dag S^{-1}s.

    A pseudoinverse is used only when the existing basis is already numerically
    singular at the requested eigenvalue floor.
    """
    if len(basis)==0:
        return 1.0

    S=nuclear_overlap_matrix(basis) if overlap_matrix is None else np.asarray(
        overlap_matrix,dtype=complex
    )
    s=np.array([
        gaussian_overlap_general(
            b.q,b.p,b.A,
            candidate.q,candidate.p,candidate.A,
        )
        for b in basis
    ],dtype=complex)

    eig=np.linalg.eigvalsh(S).real
    if np.min(eig)>float(eigenvalue_floor):
        coeff=np.linalg.solve(S,s)
    else:
        coeff=np.linalg.pinv(
            S,
            rcond=float(eigenvalue_floor)/max(np.max(eig),1.0),
        )@s

    value=1.0-np.real(np.vdot(s,coeff))
    return float(max(value,0.0))


def residual_capture_gain(
    candidate,
    residual,
    points,
    dx,
    basis,
    overlap_matrix=None,
    orthogonal_norm_floor=1e-8,
):
    r"""One-step Hilbert residual reduction available from one Gaussian pair.

    The spinor-complete candidate subspace is

        span{g_c |d_0>, g_c |d_1>}.

    Let
        g_c^\perp = (I-P_B)g_c
    be the component orthogonal to the current *nuclear* Gaussian span.  Because the
    current spinor-complete projection residual is orthogonal to every existing
    g_i|d_a>,

        <g_c^\perp d_a|r> = <g_c d_a|r>.

    The exact squared residual reduction obtained by adding the candidate electronic
    pair and reprojecting is therefore

        Delta = sum_a |<g_c|r_a>|^2 / ||g_c^\perp||^2.

    This is the v0.13 residual-greedy score.
    """
    residual=np.asarray(residual,dtype=complex)
    points=np.asarray(points,float)
    if residual.shape!=points.shape[:-1]+(2,):
        raise ValueError("residual and points have incompatible shapes.")

    nperp=candidate_orthogonal_norm(
        candidate,
        basis,
        overlap_matrix=overlap_matrix,
    )
    if nperp<float(orthogonal_norm_floor):
        return ResidualCandidateScore(
            candidate_index=-1,
            predicted_gain=0.0,
            orthogonal_norm=nperp,
        )

    g=gaussian_nd(
        points,
        candidate.q,
        candidate.p,
        candidate.A,
    )
    b=np.array([
        np.vdot(g,residual[...,a])*float(dx)*float(dx)
        for a in range(2)
    ])
    gain=float(np.sum(np.abs(b)**2)/nperp)

    return ResidualCandidateScore(
        candidate_index=-1,
        predicted_gain=max(gain,0.0),
        orthogonal_norm=nperp,
    )


def _candidate_key(candidate):
    return (
        tuple(np.round(candidate.q,12)),
        tuple(np.round(candidate.p,12)),
        tuple(np.round(candidate.A.reshape(-1),12)),
        int(candidate.state),
    )


def _basis_keys(basis):
    return {
        (
            tuple(np.round(np.asarray(b.q,float),12)),
            tuple(np.round(np.asarray(b.p,float),12)),
            tuple(np.round(np.asarray(b.A,float).reshape(-1),12)),
            int(b.state),
        )
        for b in basis
    }


def rank_residual_candidates(
    residual,
    points,
    dx,
    basis,
    candidates,
    condition_limit=1e8,
    orthogonal_norm_floor=1e-8,
):
    """Rank candidates by exact one-step spinor-complete residual capture."""
    S=nuclear_overlap_matrix(basis)
    keys=_basis_keys(basis)
    ranked=[]

    for idx,candidate in enumerate(candidates):
        if _candidate_key(candidate) in keys:
            continue

        score=residual_capture_gain(
            candidate,
            residual,
            points,
            dx,
            basis,
            overlap_matrix=S,
            orthogonal_norm_floor=orthogonal_norm_floor,
        )
        if score.orthogonal_norm<orthogonal_norm_floor:
            continue

        # Conditioning check on the proposed expanded nuclear basis.
        trial=basis+[candidate.to_tbf(-1)]
        cond=float(np.linalg.cond(nuclear_overlap_matrix(trial)))
        if not np.isfinite(cond) or cond>float(condition_limit):
            continue

        ranked.append(
            ResidualCandidateScore(
                candidate_index=idx,
                predicted_gain=score.predicted_gain,
                orthogonal_norm=score.orthogonal_norm,
            )
        )

    ranked.sort(
        key=lambda item:(-item.predicted_gain,item.candidate_index)
    )
    return ranked


def build_residual_greedy_basis(
    target_psi,
    points,
    dx,
    provider,
    seed_basis,
    candidates,
    max_basis=11,
    top_k_density_screen=1,
    condition_limit=1e5,
    orthogonal_norm_floor=1e-8,
    minimum_gain=1e-10,
    density_screen=True,
):
    r"""Build an initial Gaussian bank by deterministic residual reduction.

    Pure residual-greedy mode
    -------------------------
    `density_screen=False` or `top_k_density_screen=1` chooses the candidate with the
    largest exact one-step Hilbert-space residual gain.

    Observable-aware screened mode
    --------------------------------
    When `density_screen=True` and `top_k_density_screen>1`, the algorithm first keeps
    only the top-K candidates according to the rigorous Hilbert residual gain.  It
    then trial-projects those candidates and selects the one with the smallest
    *initial reduced electronic density* error.

    The density criterion is therefore a second-stage screen, not a replacement for
    residual reduction.  Ties are resolved by smaller wavefunction residual, then
    smaller condition number, then larger predicted gain.
    """
    target_psi=np.asarray(target_psi,dtype=complex)
    points=np.asarray(points,float)

    basis=[
        DynamicGraphTBF(
            uid=int(b.uid),
            state=int(b.state),
            q=np.asarray(b.q,float).copy(),
            p=np.asarray(b.p,float).copy(),
            A=np.asarray(b.A,float).copy(),
            node=b.node,
            spawned_targets=set(getattr(b,"spawned_targets",set())),
        )
        for b in seed_basis
    ]
    if not basis:
        raise ValueError("seed_basis must contain at least one Gaussian.")
    if int(max_basis)<len(basis):
        raise ValueError("max_basis cannot be smaller than seed basis.")

    target_density=normalized_grid_density(target_psi,dx)
    history=[]
    chosen=[]

    projection=project_grid_wavefunction_to_spinor_complete_basis(
        target_psi,
        points,
        dx,
        basis,
        provider,
    )

    while len(basis)<int(max_basis):
        residual=target_psi-projection.projected_wavefunction
        ranked=rank_residual_candidates(
            residual,
            points,
            dx,
            basis,
            candidates,
            condition_limit=condition_limit,
            orthogonal_norm_floor=orthogonal_norm_floor,
        )

        if not ranked:
            break
        if ranked[0].predicted_gain<float(minimum_gain):
            break

        k=max(1,min(int(top_k_density_screen),len(ranked)))
        screened=ranked[:k]

        if density_screen and k>1:
            trial_rows=[]
            for score in screened:
                candidate=candidates[score.candidate_index]
                trial_basis=basis+[
                    candidate.to_tbf(len(basis),"residual_trial")
                ]
                trial_projection=project_grid_wavefunction_to_spinor_complete_basis(
                    target_psi,
                    points,
                    dx,
                    trial_basis,
                    provider,
                )
                rho=normalized_grid_density(
                    trial_projection.projected_wavefunction,
                    dx,
                )
                density_error=float(np.linalg.norm(
                    rho-target_density,
                    ord="fro",
                ))
                trial_rows.append((
                    density_error,
                    float(trial_projection.relative_residual),
                    float(trial_projection.condition_number),
                    -float(score.predicted_gain),
                    int(score.candidate_index),
                    score,
                    trial_projection,
                ))

            _,_,_,_,_,selected_score,new_projection=min(trial_rows)
        else:
            selected_score=screened[0]
            candidate=candidates[selected_score.candidate_index]
            trial_basis=basis+[
                candidate.to_tbf(len(basis),"residual_trial")
            ]
            new_projection=project_grid_wavefunction_to_spinor_complete_basis(
                target_psi,
                points,
                dx,
                trial_basis,
                provider,
            )

        candidate=candidates[selected_score.candidate_index]
        old_residual=float(projection.relative_residual)

        basis.append(
            candidate.to_tbf(len(basis),"residual_selected")
        )
        chosen.append(int(selected_score.candidate_index))
        projection=new_projection

        rho=normalized_grid_density(
            projection.projected_wavefunction,
            dx,
        )
        density_error=float(np.linalg.norm(
            rho-target_density,
            ord="fro",
        ))
        actual_reduction=old_residual-float(projection.relative_residual)

        history.append(
            ResidualSelectionStep(
                basis_size=len(basis),
                selected_candidate_index=int(selected_score.candidate_index),
                selected_label=candidate.label,
                predicted_gain=float(selected_score.predicted_gain),
                actual_residual_reduction=float(actual_reduction),
                projection_fidelity=float(projection.fidelity),
                relative_residual=float(projection.relative_residual),
                density_error=density_error,
                condition_number=float(projection.condition_number),
            )
        )

    return ResidualBasisBuild(
        basis=basis,
        projection=projection,
        history=history,
        target_density=target_density,
        candidate_indices=chosen,
    )


@dataclass(frozen=True)
class PreparedGaussianDictionary:
    candidates: list
    grid_values: np.ndarray
    grid_norms: np.ndarray
    dx: float
    grid_shape: tuple


def prepare_gaussian_dictionary(
    candidates,
    points,
    dx,
):
    """Precompute candidate Gaussian values for fast residual screening."""
    points=np.asarray(points,float)
    shape=points.shape[:-1]
    G=np.empty((len(candidates),int(np.prod(shape))),dtype=complex)

    for i,candidate in enumerate(candidates):
        G[i]=gaussian_nd(
            points,
            candidate.q,candidate.p,candidate.A,
        ).reshape(-1)

    area=float(dx)*float(dx)
    norms=np.real(np.sum(np.abs(G)**2,axis=1)*area)

    return PreparedGaussianDictionary(
        candidates=list(candidates),
        grid_values=G,
        grid_norms=norms,
        dx=float(dx),
        grid_shape=tuple(shape),
    )


def _grid_basis_values(basis,points):
    points=np.asarray(points,float)
    return np.asarray([
        gaussian_nd(
            points,b.q,b.p,b.A
        ).reshape(-1)
        for b in basis
    ],dtype=complex)


def build_residual_greedy_basis_prepared(
    target_psi,
    points,
    dx,
    provider,
    seed_basis,
    prepared_dictionary,
    max_basis=11,
    top_k_density_screen=30,
    condition_limit=1e5,
    orthogonal_norm_floor=1e-8,
    minimum_gain=1e-10,
    density_screen=True,
):
    """Vectorized version of `build_residual_greedy_basis`.

    The mathematics is identical, but all candidate wavefunctions are precomputed and
    the one-step Galerkin residual gains are evaluated simultaneously with dense NumPy
    linear algebra.  Only the selected candidate is reprojected with the full v0.12
    projection routine.
    """
    target_psi=np.asarray(target_psi,dtype=complex)
    points=np.asarray(points,float)
    prepared=prepared_dictionary

    if tuple(points.shape[:-1])!=prepared.grid_shape:
        raise ValueError("prepared dictionary grid shape mismatch.")
    if abs(float(dx)-prepared.dx)>1e-15:
        raise ValueError("prepared dictionary dx mismatch.")

    basis=[
        DynamicGraphTBF(
            uid=int(b.uid),
            state=int(b.state),
            q=np.asarray(b.q,float).copy(),
            p=np.asarray(b.p,float).copy(),
            A=np.asarray(b.A,float).copy(),
            node=b.node,
            spawned_targets=set(getattr(b,"spawned_targets",set())),
        )
        for b in seed_basis
    ]
    if not basis:
        raise ValueError("seed_basis cannot be empty.")

    target_density=normalized_grid_density(target_psi,dx)
    target_norm=float(
        np.sum(np.abs(target_psi)**2)
        *float(dx)*float(dx)
    )
    area=float(dx)*float(dx)
    Rshape=target_psi.shape
    history=[]
    chosen=[]

    projection=project_grid_wavefunction_to_spinor_complete_basis(
        target_psi,points,dx,basis,provider
    )

    Gcand=prepared.grid_values
    candidates=prepared.candidates

    while len(basis)<int(max_basis):
        Gbasis=_grid_basis_values(basis,points)
        S=(Gbasis.conj()@Gbasis.T)*area
        S=0.5*(S+S.conj().T)
        Sinv=np.linalg.pinv(S,rcond=1e-12)

        residual=(
            target_psi-projection.projected_wavefunction
        ).reshape(-1,2)

        # s[:,c] = <g_i|g_c>
        s=(Gbasis.conj()@Gcand.T)*area
        alpha=Sinv@s

        nperp=prepared.grid_norms-np.real(
            np.sum(np.conj(s)*alpha,axis=0)
        )
        nperp=np.maximum(nperp,0.0)

        b=(Gcand.conj()@residual)*area
        rb=(Gbasis.conj()@residual)*area

        # <g_c^perp|r> = <g_c|r> - alpha^dag <g_basis|r>
        bperp=b-np.einsum(
            "ic,ia->ca",
            np.conj(alpha),
            rb,
        )

        gains=np.full(len(candidates),-np.inf)
        good=nperp>=float(orthogonal_norm_floor)
        gains[good]=(
            np.sum(np.abs(bperp[good])**2,axis=1)
            /nperp[good]
        )

        # Block exact duplicates/current-span elements.
        current_keys=_basis_keys(basis)
        for idx,candidate in enumerate(candidates):
            if _candidate_key(candidate) in current_keys:
                gains[idx]=-np.inf

        finite=np.flatnonzero(np.isfinite(gains))
        if len(finite)==0:
            break

        order=finite[np.argsort(gains[finite])[::-1]]
        if gains[order[0]]<float(minimum_gain):
            break

        k=max(1,min(int(top_k_density_screen),len(order)))
        top=order[:k]
        trial_rows=[]

        for idx in top:
            a=alpha[:,idx]
            gperp=Gcand[idx]-a@Gbasis
            coeff=bperp[idx]/nperp[idx]

            trial_flat=(
                projection.projected_wavefunction.reshape(-1,2)
                +gperp[:,None]*coeff[None,:]
            )
            trial_psi=trial_flat.reshape(Rshape)
            rho=normalized_grid_density(trial_psi,dx)
            density_error=float(np.linalg.norm(
                rho-target_density,ord="fro"
            ))
            rel_residual=float(
                projection.relative_residual
                -gains[idx]/max(target_norm,1e-30)
            )

            expanded=np.empty(
                (len(basis)+1,len(basis)+1),
                dtype=complex,
            )
            expanded[:-1,:-1]=S
            expanded[:-1,-1]=s[:,idx]
            expanded[-1,:-1]=np.conj(s[:,idx])
            expanded[-1,-1]=prepared.grid_norms[idx]
            cond=float(np.linalg.cond(expanded))
            if not np.isfinite(cond) or cond>condition_limit:
                continue

            if density_screen and k>1:
                key=(
                    density_error,
                    rel_residual,
                    cond,
                    -float(gains[idx]),
                    int(idx),
                )
            else:
                key=(
                    -float(gains[idx]),
                    density_error,
                    rel_residual,
                    cond,
                    int(idx),
                )
            trial_rows.append((key,int(idx)))

        if not trial_rows:
            break

        _,idx=min(trial_rows,key=lambda row:row[0])
        candidate=candidates[idx]
        old_residual=float(projection.relative_residual)

        basis.append(
            candidate.to_tbf(
                len(basis),
                "residual_selected_fast",
            )
        )
        chosen.append(int(idx))

        projection=project_grid_wavefunction_to_spinor_complete_basis(
            target_psi,points,dx,basis,provider
        )
        rho=normalized_grid_density(
            projection.projected_wavefunction,dx
        )
        density_error=float(np.linalg.norm(
            rho-target_density,ord="fro"
        ))

        history.append(
            ResidualSelectionStep(
                basis_size=len(basis),
                selected_candidate_index=int(idx),
                selected_label=candidate.label,
                predicted_gain=float(gains[idx]),
                actual_residual_reduction=float(
                    old_residual-projection.relative_residual
                ),
                projection_fidelity=float(projection.fidelity),
                relative_residual=float(projection.relative_residual),
                density_error=density_error,
                condition_number=float(projection.condition_number),
            )
        )

    return ResidualBasisBuild(
        basis=basis,
        projection=projection,
        history=history,
        target_density=target_density,
        candidate_indices=chosen,
    )
