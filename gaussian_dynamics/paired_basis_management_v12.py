from dataclasses import dataclass
import numpy as np

from .basis_management import overlap_conditioning


@dataclass
class PairedPruningResult:
    keep: np.ndarray
    removed: tuple
    coefficients_matrix: np.ndarray
    nuclear_overlap: np.ndarray
    projection_loss: float
    condition_before: float
    condition_after: float


def spinor_wavefunction_norm(Cmat, Snuc):
    C=np.asarray(Cmat,dtype=complex)
    S=np.asarray(Snuc,dtype=complex)
    if C.ndim!=2 or S.shape!=(len(C),len(C)):
        raise ValueError("incompatible Cmat/Snuc shapes.")
    return float(sum(
        np.real(np.vdot(C[:,a],S@C[:,a]))
        for a in range(C.shape[1])
    ))


def project_spinor_coefficients_to_subset(Cmat, Snuc, keep):
    """Project every electronic component onto one retained nuclear Gaussian subset."""
    C=np.asarray(Cmat,dtype=complex)
    S=np.asarray(Snuc,dtype=complex)
    keep=np.asarray(keep,dtype=int)

    if C.ndim!=2 or S.shape!=(len(C),len(C)):
        raise ValueError("incompatible Cmat/Snuc shapes.")
    if len(keep)==0:
        raise ValueError("cannot project onto an empty nuclear basis.")

    Sk=S[np.ix_(keep,keep)]
    Sall=S[np.ix_(keep,np.arange(len(C)))]

    Cnew=np.zeros((len(keep),C.shape[1]),dtype=complex)
    for a in range(C.shape[1]):
        Cnew[:,a]=np.linalg.solve(
            Sk,
            Sall@C[:,a],
        )

    old=spinor_wavefunction_norm(C,S)
    new=spinor_wavefunction_norm(Cnew,Sk)
    return Cnew,max(old-new,0.0)


def _candidate_smallest_mode(S):
    eig,U=np.linalg.eigh(np.asarray(S,dtype=complex))
    vec=U[:,np.argmin(eig.real)]
    return int(np.argmax(np.abs(vec)))


def prune_nuclear_gaussian_pairs(
    Cmat,
    Snuc,
    condition_limit=1e9,
    eigenvalue_floor=1e-10,
    max_projection_loss=1e-7,
):
    """Prune whole nuclear Gaussians while retaining both electronic components."""
    C=np.asarray(Cmat,dtype=complex).copy()
    S=np.asarray(Snuc,dtype=complex).copy()
    keep_global=np.arange(len(C),dtype=int)
    removed=[]
    total_loss=0.0

    before=overlap_conditioning(S).condition_number

    while len(C)>1:
        report=overlap_conditioning(S)
        need=(
            report.condition_number>condition_limit
            or report.smallest_eigenvalue<eigenvalue_floor
        )
        if not need:
            break

        candidate=_candidate_smallest_mode(S)
        local_keep=np.array(
            [i for i in range(len(C)) if i!=candidate],
            dtype=int,
        )
        Ctrial,loss=project_spinor_coefficients_to_subset(
            C,S,local_keep
        )
        if total_loss+loss>max_projection_loss:
            break

        removed.append(int(keep_global[candidate]))
        total_loss+=loss
        keep_global=keep_global[local_keep]
        C=Ctrial
        S=S[np.ix_(local_keep,local_keep)]

    after=overlap_conditioning(S).condition_number

    return PairedPruningResult(
        keep=keep_global,
        removed=tuple(removed),
        coefficients_matrix=C,
        nuclear_overlap=S,
        projection_loss=float(total_loss),
        condition_before=float(before),
        condition_after=float(after),
    )
