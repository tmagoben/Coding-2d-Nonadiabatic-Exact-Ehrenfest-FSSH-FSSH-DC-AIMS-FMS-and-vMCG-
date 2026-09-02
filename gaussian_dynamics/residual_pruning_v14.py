from dataclasses import dataclass
import numpy as np

from .paired_basis_management_v12 import (
    project_spinor_coefficients_to_subset,
    spinor_wavefunction_norm,
)


@dataclass(frozen=True)
class LeaveOneOutScore:
    index: int
    uid: int
    orthogonal_norm: float
    absolute_projection_loss: float
    fractional_projection_loss: float


@dataclass
class LowLossPruningResult:
    keep: np.ndarray
    removed_index: int
    removed_uid: int
    coefficients_matrix: np.ndarray
    nuclear_overlap: np.ndarray
    absolute_projection_loss: float
    fractional_projection_loss: float
    condition_before: float
    condition_after: float


def leave_one_out_projection_losses(Cmat, Snuc, uids=None):
    r"""Exact one-Gaussian deletion loss for the represented spinor wavefunction.

    Let `S` be the nuclear Gaussian Gram matrix and let `C_j` be the complete
    electronic coefficient row carried by Gaussian j.  The component of g_j
    orthogonal to the span of all remaining Gaussians has norm

        n_j = 1 / (S^{-1})_jj.

    Deleting Gaussian j and optimally projecting the old wavefunction into the
    remaining basis loses exactly

        L_j = n_j * sum_a |C_ja|^2.

    One dense solve/inverse supplies every `n_j`, so all N deletion scores cost
    O(N^3 + Ns), rather than N separate O(N^3) projections.
    """
    C=np.asarray(Cmat,dtype=complex)
    S=np.asarray(Snuc,dtype=complex)

    if C.ndim!=2 or S.shape!=(len(C),len(C)):
        raise ValueError("incompatible Cmat/Snuc shapes.")
    if len(C)<2:
        return []

    if uids is None:
        uids=np.arange(len(C))
    uids=np.asarray(uids)
    if uids.shape!=(len(C),):
        raise ValueError("uids must have one entry per Gaussian.")

    Sinv=np.linalg.solve(S,np.eye(len(S),dtype=complex))
    diag=np.real(np.diag(Sinv))

    old_norm=spinor_wavefunction_norm(C,S)
    scores=[]

    for j in range(len(C)):
        if diag[j]<=0.0:
            continue
        nperp=float(1.0/diag[j])
        row_weight=float(np.sum(np.abs(C[j])**2))
        loss=float(max(nperp*row_weight,0.0))
        frac=float(loss/max(old_norm,1e-30))
        scores.append(
            LeaveOneOutScore(
                index=int(j),
                uid=int(uids[j]),
                orthogonal_norm=nperp,
                absolute_projection_loss=loss,
                fractional_projection_loss=frac,
            )
        )

    scores.sort(
        key=lambda x:(x.fractional_projection_loss,x.index)
    )
    return scores


def prune_low_loss_gaussian_pair(
    Cmat,
    Snuc,
    uids,
    max_fractional_loss=1e-7,
    protected_uids=(),
    require_condition_improvement=False,
):
    """Remove the lowest-loss whole nuclear Gaussian/electronic pair, if admissible."""
    C=np.asarray(Cmat,dtype=complex)
    S=np.asarray(Snuc,dtype=complex)
    uids=np.asarray(uids)
    protected={int(x) for x in protected_uids}

    if len(C)<=1:
        return None

    before=float(np.linalg.cond(S))
    scores=leave_one_out_projection_losses(
        C,S,uids=uids
    )

    for score in scores:
        if score.uid in protected:
            continue
        if score.fractional_projection_loss>float(max_fractional_loss):
            return None

        keep=np.array(
            [i for i in range(len(C)) if i!=score.index],
            dtype=int,
        )
        Snew=S[np.ix_(keep,keep)]
        after=float(np.linalg.cond(Snew))

        if require_condition_improvement and not after<before:
            continue

        Cnew,actual_loss=project_spinor_coefficients_to_subset(
            C,S,keep
        )
        old_norm=spinor_wavefunction_norm(C,S)
        actual_frac=float(actual_loss/max(old_norm,1e-30))

        return LowLossPruningResult(
            keep=keep,
            removed_index=int(score.index),
            removed_uid=int(score.uid),
            coefficients_matrix=Cnew,
            nuclear_overlap=Snew,
            absolute_projection_loss=float(actual_loss),
            fractional_projection_loss=actual_frac,
            condition_before=before,
            condition_after=after,
        )

    return None
