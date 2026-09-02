from collections import Counter
import numpy as np


def generation_histogram(lineage, active_basis=None):
    if active_basis is None:
        uids=set(lineage)
    else:
        uids={int(b.uid) for b in active_basis}

    counts=Counter(
        int(info["generation"])
        for uid,info in lineage.items()
        if int(uid) in uids
    )
    return dict(sorted(counts.items()))


def width_determinants(basis):
    return np.asarray([np.linalg.det(np.asarray(b.A,float)) for b in basis],float)


def width_diversity_ratio(basis):
    dets=width_determinants(basis)
    if len(dets)==0:
        return 1.0
    return float(np.max(dets)/np.min(dets))


def overlap_spectrum_metrics(S, relative_cutoff=1e-12):
    S=np.asarray(S,dtype=complex)
    eig=np.linalg.eigvalsh(0.5*(S+S.conj().T)).real
    maxeig=max(float(np.max(eig)),0.0)
    cutoff=relative_cutoff*max(maxeig,1.0)
    retained=eig[eig>cutoff]

    rank=len(retained)
    cond=np.inf if rank<len(eig) else float(np.max(retained)/np.min(retained))

    if rank==0:
        entropy=0.0
        spectral_effective_rank=0.0
    else:
        p=retained/np.sum(retained)
        entropy=float(-np.sum(p*np.log(p)))
        spectral_effective_rank=float(np.exp(entropy))

    return {
        "eigenvalues":eig,
        "numerical_rank":int(rank),
        "condition_number":float(cond),
        "spectral_entropy":entropy,
        "spectral_effective_rank":spectral_effective_rank,
    }


def canonical_coefficient_weights(C,S,relative_cutoff=1e-12):
    """Weights of Psi in the canonical orthonormalized overlap eigenbasis."""
    C=np.asarray(C,dtype=complex)
    S=np.asarray(S,dtype=complex)
    eig,U=np.linalg.eigh(0.5*(S+S.conj().T))

    cutoff=relative_cutoff*max(float(np.max(eig.real)),1.0)
    mask=eig.real>cutoff
    eig_r=eig[mask].real
    U_r=U[:,mask]

    # d = sqrt(s) U^dag C because chi = Phi U s^-1/2 and Psi=Phi C.
    d=np.sqrt(eig_r)*(U_r.conj().T@C)
    weights=np.abs(d)**2
    norm=float(np.sum(weights))

    if norm>0:
        probabilities=weights/norm
        participation=1.0/float(np.sum(probabilities**2))
    else:
        probabilities=np.zeros_like(weights)
        participation=0.0

    return {
        "canonical_amplitudes":d,
        "weights":weights,
        "probabilities":probabilities,
        "norm":norm,
        "participation_ratio":participation,
        "retained_rank":int(len(eig_r)),
    }


def lineage_depth(lineage):
    return max(
        (int(info["generation"]) for info in lineage.values()),
        default=0,
    )


def basis_completeness_report(run):
    S=np.asarray(run["final_overlap"],dtype=complex)
    basis=run["final_basis"]
    lineage=run["lineage"]

    spectrum=overlap_spectrum_metrics(S)
    coeff=canonical_coefficient_weights(
        run["final_coefficients"],
        S,
    )

    return {
        "basis_size":len(basis),
        "lineage_depth":lineage_depth(lineage),
        "generation_histogram":generation_histogram(lineage,basis),
        "width_determinants":width_determinants(basis),
        "width_diversity_ratio":width_diversity_ratio(basis),
        "overlap_spectrum":spectrum,
        "canonical_participation_ratio":coeff["participation_ratio"],
        "canonical_probabilities":coeff["probabilities"],
    }
