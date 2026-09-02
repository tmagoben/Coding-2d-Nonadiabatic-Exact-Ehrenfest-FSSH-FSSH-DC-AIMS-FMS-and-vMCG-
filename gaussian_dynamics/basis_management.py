from dataclasses import dataclass
import numpy as np


@dataclass
class BasisConditioningReport:
    eigenvalues: np.ndarray
    condition_number: float
    numerical_rank: int
    smallest_eigenvalue: float


@dataclass
class PruningResult:
    keep: np.ndarray
    removed: tuple
    coefficients: np.ndarray
    overlap: np.ndarray
    projection_loss: float
    condition_before: float
    condition_after: float


def overlap_conditioning(S, relative_floor=1e-12):
    S = np.asarray(S, dtype=complex)
    if S.ndim != 2 or S.shape[0] != S.shape[1]:
        raise ValueError("S must be square")
    if not np.allclose(S, S.conj().T, atol=1e-10):
        raise ValueError("S must be Hermitian")

    eig = np.linalg.eigvalsh(S).real
    maxeig = max(float(np.max(eig)), 0.0)
    floor = relative_floor * max(maxeig, 1.0)
    positive = eig[eig > floor]
    rank = len(positive)

    if rank == 0:
        cond = np.inf
    elif rank < len(eig):
        cond = np.inf
    else:
        cond = float(np.max(eig) / np.min(eig))

    return BasisConditioningReport(
        eigenvalues=eig,
        condition_number=cond,
        numerical_rank=rank,
        smallest_eigenvalue=float(np.min(eig)),
    )


def project_coefficients_to_subset(C, S, keep, regularization=0.0):
    """Least-squares Hilbert-space projection onto a retained basis subset.

    For Psi_old = sum_j C_j phi_j and retained indices K, solve

        S_KK C_new = S_K,all C_old.

    The returned projection_loss is ||Psi_old-P_K Psi_old||^2.
    """
    C = np.asarray(C, dtype=complex)
    S = np.asarray(S, dtype=complex)
    keep = np.asarray(keep, dtype=int)

    if S.shape != (len(C), len(C)):
        raise ValueError("C and S shapes are incompatible")
    if len(keep) == 0:
        raise ValueError("cannot project onto an empty basis")

    Sk = S[np.ix_(keep, keep)].copy()
    if regularization > 0.0:
        Sk += regularization * np.eye(len(keep))
    b = S[np.ix_(keep, np.arange(len(C)))] @ C
    Cnew = np.linalg.solve(Sk, b)

    old_norm = float(np.real(np.vdot(C, S @ C)))
    projected_norm = float(np.real(np.vdot(Cnew, S[np.ix_(keep, keep)] @ Cnew)))
    loss = max(old_norm - projected_norm, 0.0)

    return Cnew, loss


def _candidate_from_smallest_mode(S, protected):
    eig, U = np.linalg.eigh(np.asarray(S, dtype=complex))
    vec = U[:, np.argmin(eig.real)]
    order = np.argsort(np.abs(vec))[::-1]
    protected = set(int(i) for i in protected)
    for idx in order:
        if int(idx) not in protected:
            return int(idx)
    return None


def prune_redundant_basis(
    C,
    S,
    condition_limit=1e8,
    eigenvalue_floor=1e-9,
    max_projection_loss=1e-8,
    protected_indices=(),
):
    """Iteratively remove near-redundant basis functions with controlled wavefunction loss."""
    Cwork = np.asarray(C, dtype=complex).copy()
    Swork = np.asarray(S, dtype=complex).copy()
    keep_global = np.arange(len(Cwork), dtype=int)
    protected_global = set(int(i) for i in protected_indices)
    removed = []
    total_loss = 0.0

    before = overlap_conditioning(Swork).condition_number

    while len(Cwork) > 1:
        report = overlap_conditioning(Swork)
        need_prune = (
            report.condition_number > condition_limit
            or report.smallest_eigenvalue < eigenvalue_floor
        )
        if not need_prune:
            break

        protected_local = [
            local for local, global_idx in enumerate(keep_global)
            if int(global_idx) in protected_global
        ]
        candidate = _candidate_from_smallest_mode(Swork, protected_local)
        if candidate is None:
            break

        local_keep = np.array([i for i in range(len(Cwork)) if i != candidate], dtype=int)
        Ctrial, loss = project_coefficients_to_subset(Cwork, Swork, local_keep)

        if total_loss + loss > max_projection_loss:
            break

        removed.append(int(keep_global[candidate]))
        total_loss += loss
        keep_global = keep_global[local_keep]
        Swork = Swork[np.ix_(local_keep, local_keep)]
        Cwork = Ctrial

    after = overlap_conditioning(Swork).condition_number

    return PruningResult(
        keep=keep_global,
        removed=tuple(removed),
        coefficients=Cwork,
        overlap=Swork,
        projection_loss=float(total_loss),
        condition_before=float(before),
        condition_after=float(after),
    )


def canonical_orthogonalizer(S, relative_cutoff=1e-10):
    """Return X with X^dagger S X = I in the retained canonical subspace."""
    S = np.asarray(S, dtype=complex)
    eig, U = np.linalg.eigh(S)
    cutoff = relative_cutoff * max(float(np.max(eig.real)), 1.0)
    mask = eig.real > cutoff
    if not np.any(mask):
        raise np.linalg.LinAlgError("overlap matrix has no retained positive eigenspace")
    X = U[:, mask] @ np.diag(1.0 / np.sqrt(eig[mask].real))
    return X, eig, mask
