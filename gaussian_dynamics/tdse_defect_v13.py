from dataclasses import dataclass
import numpy as np

from .gaussian_nd import gaussian_nd, gaussian_nd_time_derivative
from .spinor_complete_lvc_v12 import (
    build_spinor_complete_lvc_matrices,
    build_spinor_complete_time_matrix,
    coefficients_matrix,
)
from .dynamic_graph_aims import DynamicGraphTBF, _kinematics
from .born_huang_grid_v12 import apply_spectral_kinetic_to_basis_fields
from .residual_basis_v13 import (
    nuclear_overlap_matrix,
    GaussianCandidate,
)


@dataclass(frozen=True)
class TDSEDefect:
    residual: np.ndarray
    wavefunction: np.ndarray
    wavefunction_time_derivative: np.ndarray
    hamiltonian_wavefunction: np.ndarray
    residual_norm: float
    relative_to_hpsi: float
    coefficient_derivative: np.ndarray
    projected_residual_norm: float


@dataclass(frozen=True)
class DefectCandidateScore:
    candidate_index: int
    captured_defect_norm: float
    capture_fraction: float
    orthogonal_norm: float
    expanded_condition_number: float


def reconstruct_spinor_complete_wavefunction(
    flat_coefficients,
    basis,
    points,
):
    C=coefficients_matrix(flat_coefficients,len(basis))
    points=np.asarray(points,float)
    psi=np.zeros(points.shape[:-1]+(2,),dtype=complex)

    for i,b in enumerate(basis):
        g=gaussian_nd(points,b.q,b.p,b.A)
        psi += g[...,None]*C[i][None,None,:]

    return psi


def spinor_complete_coefficient_derivative(
    flat_coefficients,
    basis,
    provider,
):
    """Continuous Galerkin coefficient derivative for the moving Gaussian basis."""
    C=np.asarray(flat_coefficients,dtype=complex)
    S,H,_=build_spinor_complete_lvc_matrices(
        basis,provider
    )

    qdots=[]
    pdots=[]
    for b in basis:
        qdot,pdot=_kinematics(b,provider)
        qdots.append(qdot)
        pdots.append(pdot)

    qdots=np.asarray(qdots,float)
    pdots=np.asarray(pdots,float)

    T=build_spinor_complete_time_matrix(
        basis,qdots,pdots
    )
    Cdot=np.linalg.solve(
        S,
        -(1j*H+T)@C,
    )
    return Cdot,qdots,pdots,S,H,T


def reconstruct_spinor_complete_time_derivative(
    flat_coefficients,
    coefficient_derivative,
    basis,
    points,
    qdots,
    pdots,
):
    C=coefficients_matrix(flat_coefficients,len(basis))
    Cdot=coefficients_matrix(
        coefficient_derivative,len(basis)
    )
    points=np.asarray(points,float)

    psidot=np.zeros(
        points.shape[:-1]+(2,),
        dtype=complex,
    )

    for i,b in enumerate(basis):
        g=gaussian_nd(points,b.q,b.p,b.A)
        gdot=gaussian_nd_time_derivative(
            points,
            b.q,b.p,b.A,
            qdots[i],pdots[i],
        )
        psidot += (
            gdot[...,None]*C[i][None,None,:]
            +g[...,None]*Cdot[i][None,None,:]
        )

    return psidot


def apply_lvc_grid_hamiltonian(psi,grid):
    """Apply the same periodic FFT kinetic + global diabatic LVC potential."""
    psi=np.asarray(psi,dtype=complex)
    if psi.shape!=grid.points.shape[:-1]+(2,):
        raise ValueError("psi has incompatible shape for this grid.")

    Tpsi=apply_spectral_kinetic_to_basis_fields(
        psi[None,...],
        grid,
    )[0]
    Vpsi=np.einsum(
        "...ab,...b->...a",
        grid.diabatic_potential,
        psi,
    )
    return Tpsi+Vpsi


def _grid_inner(a,b,dx):
    return np.vdot(
        np.asarray(a,dtype=complex).reshape(-1),
        np.asarray(b,dtype=complex).reshape(-1),
    )*float(dx)*float(dx)


def compute_tdse_defect(
    flat_coefficients,
    basis,
    provider,
    grid,
):
    r"""Compute the instantaneous Schrödinger/Galerkin defect.

    The residual is

        R = i dPsi/dt - H Psi.

    For an exact solution R=0.  For a finite moving Gaussian Galerkin basis, the
    projected equations make R orthogonal to the represented tangent directions, but
    components outside the finite span remain.

    v0.13 uses the norm of that component as a basis-incompleteness diagnostic.
    """
    Cdot,qdots,pdots,S,H,T=spinor_complete_coefficient_derivative(
        flat_coefficients,
        basis,
        provider,
    )
    psi=reconstruct_spinor_complete_wavefunction(
        flat_coefficients,basis,grid.points
    )
    psidot=reconstruct_spinor_complete_time_derivative(
        flat_coefficients,
        Cdot,
        basis,
        grid.points,
        qdots,
        pdots,
    )
    Hpsi=apply_lvc_grid_hamiltonian(psi,grid)
    residual=1j*psidot-Hpsi

    residual_norm=float(np.sqrt(max(
        np.real(_grid_inner(residual,residual,grid.dx)),
        0.0,
    )))
    hnorm=float(np.sqrt(max(
        np.real(_grid_inner(Hpsi,Hpsi,grid.dx)),
        0.0,
    )))

    # Project residual onto the current spinor-complete basis on the diagnostic grid.
    n=len(basis)
    b=np.zeros((n,2),dtype=complex)
    for i,tbf in enumerate(basis):
        g=gaussian_nd(
            grid.points,tbf.q,tbf.p,tbf.A
        )
        for a in range(2):
            b[i,a]=np.vdot(
                g,residual[...,a]
            )*grid.area

    Snuc=nuclear_overlap_matrix(basis)
    projected=0.0
    for a in range(2):
        if np.linalg.norm(b[:,a])==0.0:
            continue
        coeff=np.linalg.lstsq(
            Snuc,b[:,a],rcond=1e-12
        )[0]
        projected += np.real(
            np.vdot(b[:,a],coeff)
        )
    projected_norm=float(np.sqrt(max(projected,0.0)))

    return TDSEDefect(
        residual=residual,
        wavefunction=psi,
        wavefunction_time_derivative=psidot,
        hamiltonian_wavefunction=Hpsi,
        residual_norm=residual_norm,
        relative_to_hpsi=float(
            residual_norm/max(hnorm,1e-30)
        ),
        coefficient_derivative=Cdot,
        projected_residual_norm=projected_norm,
    )


def defect_candidate_capture(
    candidate,
    defect,
    basis,
    grid,
    condition_limit=1e8,
    orthogonal_norm_floor=1e-8,
):
    """Residual norm captured by adding one spinor-complete Gaussian candidate."""
    if not isinstance(candidate,GaussianCandidate):
        raise TypeError("candidate must be a GaussianCandidate.")

    S=nuclear_overlap_matrix(basis)
    s=np.array([
        np.vdot(
            gaussian_nd(
                grid.points,b.q,b.p,b.A
            ),
            gaussian_nd(
                grid.points,
                candidate.q,candidate.p,candidate.A
            ),
        )*grid.area
        for b in basis
    ],dtype=complex)

    coeff=np.linalg.lstsq(S,s,rcond=1e-12)[0]
    nperp=float(max(
        1.0-np.real(np.vdot(s,coeff)),
        0.0,
    ))
    if nperp<float(orthogonal_norm_floor):
        return None

    gc=gaussian_nd(
        grid.points,
        candidate.q,candidate.p,candidate.A,
    )
    gperp=gc.copy()
    for alpha,b in zip(coeff,basis):
        gperp -= alpha*gaussian_nd(
            grid.points,b.q,b.p,b.A
        )

    bvec=np.array([
        np.vdot(gperp,defect.residual[...,a])*grid.area
        for a in range(2)
    ])
    captured=float(np.sum(np.abs(bvec)**2)/nperp)

    expanded=np.empty(
        (len(basis)+1,len(basis)+1),
        dtype=complex,
    )
    expanded[:-1,:-1]=S
    expanded[:-1,-1]=s
    expanded[-1,:-1]=np.conj(s)
    expanded[-1,-1]=1.0
    cond=float(np.linalg.cond(expanded))
    if not np.isfinite(cond) or cond>condition_limit:
        return None

    defect_sq=defect.residual_norm**2
    return DefectCandidateScore(
        candidate_index=-1,
        captured_defect_norm=float(np.sqrt(max(captured,0.0))),
        capture_fraction=float(
            captured/max(defect_sq,1e-30)
        ),
        orthogonal_norm=nperp,
        expanded_condition_number=cond,
    )


def rank_defect_candidates(
    defect,
    basis,
    candidates,
    grid,
    condition_limit=1e8,
    orthogonal_norm_floor=1e-8,
):
    ranked=[]
    for idx,candidate in enumerate(candidates):
        score=defect_candidate_capture(
            candidate,
            defect,
            basis,
            grid,
            condition_limit=condition_limit,
            orthogonal_norm_floor=orthogonal_norm_floor,
        )
        if score is None:
            continue
        ranked.append(
            DefectCandidateScore(
                candidate_index=idx,
                captured_defect_norm=score.captured_defect_norm,
                capture_fraction=score.capture_fraction,
                orthogonal_norm=score.orthogonal_norm,
                expanded_condition_number=score.expanded_condition_number,
            )
        )

    ranked.sort(
        key=lambda x:(-x.capture_fraction,x.candidate_index)
    )
    return ranked


@dataclass
class DefectEnrichmentResult:
    basis: list
    coefficients: np.ndarray
    defect_before: TDSEDefect
    defect_after: TDSEDefect
    selected_candidate_index: int
    selected_candidate: GaussianCandidate
    score: DefectCandidateScore
    actual_squared_defect_reduction: float


def enrich_basis_from_tdse_defect(
    flat_coefficients,
    basis,
    provider,
    grid,
    candidates,
    condition_limit=1e8,
    orthogonal_norm_floor=1e-8,
):
    r"""Add the candidate that captures the largest instantaneous TDSE defect.

    The new Gaussian enters with a zero two-component electronic coefficient, so the
    represented wavefunction is exactly unchanged at insertion.  What changes is the
    available tangent/Galerkin space for `Cdot`.

    For an ideal Galerkin projection, the expected decrease in squared defect is the
    candidate's captured residual norm squared.
    """
    before=compute_tdse_defect(
        flat_coefficients,basis,provider,grid
    )
    ranked=rank_defect_candidates(
        before,
        basis,
        candidates,
        grid,
        condition_limit=condition_limit,
        orthogonal_norm_floor=orthogonal_norm_floor,
    )
    if not ranked:
        return None

    score=ranked[0]
    candidate=candidates[score.candidate_index]

    new_basis=[
        DynamicGraphTBF(
            uid=int(b.uid),
            state=int(b.state),
            q=np.asarray(b.q,float).copy(),
            p=np.asarray(b.p,float).copy(),
            A=np.asarray(b.A,float).copy(),
            node=b.node,
            spawned_targets=set(getattr(b,"spawned_targets",set())),
        )
        for b in basis
    ]
    new_basis.append(
        candidate.to_tbf(
            len(new_basis),
            node_prefix="defect_selected",
        )
    )

    C=np.asarray(flat_coefficients,dtype=complex)
    new_C=np.concatenate([
        C,
        np.zeros(2,dtype=complex),
    ])

    after=compute_tdse_defect(
        new_C,new_basis,provider,grid
    )

    return DefectEnrichmentResult(
        basis=new_basis,
        coefficients=new_C,
        defect_before=before,
        defect_after=after,
        selected_candidate_index=int(score.candidate_index),
        selected_candidate=candidate,
        score=score,
        actual_squared_defect_reduction=float(
            before.residual_norm**2
            -after.residual_norm**2
        ),
    )
