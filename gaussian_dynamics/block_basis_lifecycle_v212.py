from dataclasses import dataclass
import numpy as np
from scipy import sparse


@dataclass(frozen=True)
class BlockPruneResultV212:
    basis: tuple
    coefficients: np.ndarray
    removed_index: int
    removed_uid: int
    projection_loss: float
    retained_condition_number: float


def _scalar_block_indices(n_basis, nstate, block_index):
    start=int(block_index)*int(nstate)
    return np.arange(start,start+int(nstate),dtype=int)


def insert_zero_block_v212(basis, coefficients, new_tbf, nstate, index=None):
    """Insert a Gaussian TBF with an exactly zero electronic coefficient block.

    Zero-block insertion enlarges the variational space without changing the represented
    molecular state at the instant of insertion.
    """
    basis=list(basis)
    C=np.asarray(coefficients,dtype=complex)
    s=int(nstate)
    if C.shape!=(len(basis)*s,):
        raise ValueError("coefficient dimension is inconsistent with basis*nstate.")
    if index is None:
        index=len(basis)
    index=int(index)
    if not (0<=index<=len(basis)):
        raise ValueError("insertion index is out of range.")
    basis.insert(index,new_tbf)
    pos=index*s
    Cnew=np.concatenate((C[:pos],np.zeros(s,dtype=complex),C[pos:]))
    return tuple(basis),Cnew


def prune_block_projected_v212(basis, coefficients, S, nstate, block_index):
    """Project the represented state onto the basis with one Gaussian block removed.

    Partition the metric as

        S = [[S_rr, S_rd], [S_dr, S_dd]]

    after reordering retained and deleted scalar components.  The orthogonal projection
    onto the retained span has coefficients

        C_r' = C_r + S_rr^{-1} S_rd C_d.

    The squared projection loss is

        C_d^dagger (S_dd - S_dr S_rr^{-1} S_rd) C_d.

    This is a block-level generalization of the low-loss pruning algebra and does not
    assume a spin label or a one-state-per-TBF representation.
    """
    basis=list(basis)
    C=np.asarray(coefficients,dtype=complex)
    s=int(nstate)
    n=len(basis)
    k=int(block_index)
    if not (0<=k<n):
        raise ValueError("block_index is out of range.")
    if C.shape!=(n*s,):
        raise ValueError("coefficient dimension is inconsistent with basis*nstate.")
    A=S.toarray() if sparse.issparse(S) else np.asarray(S,dtype=complex)
    if A.shape!=(n*s,n*s):
        raise ValueError("metric has incompatible shape.")

    deleted=_scalar_block_indices(n,s,k)
    mask=np.ones(n*s,dtype=bool)
    mask[deleted]=False
    retained=np.flatnonzero(mask)

    Srr=A[np.ix_(retained,retained)]
    Srd=A[np.ix_(retained,deleted)]
    Sdr=A[np.ix_(deleted,retained)]
    Sdd=A[np.ix_(deleted,deleted)]
    Cr=C[retained]
    Cd=C[deleted]

    correction=np.linalg.solve(Srr,Srd@Cd)
    Cproj=Cr+correction
    schur=Sdd-Sdr@np.linalg.solve(Srr,Srd)
    loss=float(np.real(np.vdot(Cd,schur@Cd)))
    loss=max(loss,0.0)

    kept_basis=tuple(b for i,b in enumerate(basis) if i!=k)
    return BlockPruneResultV212(
        basis=kept_basis,
        coefficients=Cproj,
        removed_index=k,
        removed_uid=int(basis[k].uid),
        projection_loss=loss,
        retained_condition_number=float(np.linalg.cond(Srr)),
    )
