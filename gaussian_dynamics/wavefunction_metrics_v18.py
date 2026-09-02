from dataclasses import dataclass
import numpy as np

from .tdse_defect_v13 import reconstruct_spinor_complete_wavefunction


def grid_inner_product(psi, phi, area):
    """Spin-summed grid inner product."""
    a=np.asarray(psi,dtype=complex)
    b=np.asarray(phi,dtype=complex)
    if a.shape!=b.shape:
        raise ValueError("wavefunctions must have equal shapes.")
    return np.vdot(a.reshape(-1),b.reshape(-1))*float(area)


def grid_wavefunction_norm(psi, area):
    value=np.real(grid_inner_product(psi,psi,area))
    return float(np.sqrt(max(value,0.0)))


def normalize_grid_wavefunction(psi, area):
    psi=np.asarray(psi,dtype=complex)
    norm=grid_wavefunction_norm(psi,area)
    if norm<=1e-15:
        raise ValueError("cannot normalize a zero grid wavefunction.")
    return psi/norm


def phase_align_wavefunction(reference, candidate, area):
    """Return candidate with the global phase minimizing its L2 distance to reference."""
    ref=normalize_grid_wavefunction(reference,area)
    cand=normalize_grid_wavefunction(candidate,area)
    z=grid_inner_product(ref,cand,area)
    if abs(z)<=1e-15:
        return cand,1.0+0.0j
    factor=np.exp(-1j*np.angle(z))
    return cand*factor,complex(factor)


def phase_aligned_fidelity(reference, candidate, area):
    ref=normalize_grid_wavefunction(reference,area)
    cand=normalize_grid_wavefunction(candidate,area)
    z=grid_inner_product(ref,cand,area)
    return float(np.clip(abs(z)**2,0.0,1.0+1e-12))


def phase_aligned_l2_error(reference, candidate, area):
    ref=normalize_grid_wavefunction(reference,area)
    aligned,_=phase_align_wavefunction(ref,candidate,area)
    diff=ref-aligned
    value=np.real(grid_inner_product(diff,diff,area))
    return float(np.sqrt(max(value,0.0)))


def nuclear_density(psi, area=None, normalize=True):
    psi=np.asarray(psi,dtype=complex)
    if psi.ndim<2:
        raise ValueError("spinor wavefunction must have a state axis.")
    rho=np.sum(np.abs(psi)**2,axis=-1)
    if normalize:
        if area is None:
            raise ValueError("area is required when normalize=True.")
        total=float(np.sum(rho)*float(area))
        if total<=1e-15:
            raise ValueError("zero nuclear-density norm.")
        rho=rho/total
    return np.asarray(rho,dtype=float)


def nuclear_density_l2_error(reference, candidate, area):
    nr=nuclear_density(reference,area,normalize=True)
    nc=nuclear_density(candidate,area,normalize=True)
    return float(np.sqrt(
        np.sum((nr-nc)**2)*float(area)
    ))


def nuclear_density_total_variation(reference, candidate, area):
    nr=nuclear_density(reference,area,normalize=True)
    nc=nuclear_density(candidate,area,normalize=True)
    return float(
        0.5*np.sum(np.abs(nr-nc))*float(area)
    )


def spatial_moments(psi, points, area):
    """Normalized nuclear centroid and covariance from a spinor grid wavefunction."""
    points=np.asarray(points,dtype=float)
    rho=nuclear_density(psi,area,normalize=True)
    if points.shape[:-1]!=rho.shape:
        raise ValueError("points and wavefunction grid shapes are incompatible.")

    coords=points.reshape(-1,points.shape[-1])
    weights=(rho.reshape(-1)*float(area))
    mean=np.sum(weights[:,None]*coords,axis=0)
    centered=coords-mean
    cov=np.einsum(
        "n,ni,nj->ij",
        weights,centered,centered,
    )
    return {
        "mean":np.asarray(mean,dtype=float),
        "covariance":np.asarray(cov,dtype=float),
        "variances":np.diag(cov).copy(),
    }


def moment_errors(reference, candidate, points, area):
    a=spatial_moments(reference,points,area)
    b=spatial_moments(candidate,points,area)
    return {
        "mean_l2":float(np.linalg.norm(
            b["mean"]-a["mean"]
        )),
        "covariance_frobenius":float(np.linalg.norm(
            b["covariance"]-a["covariance"],
            ord="fro",
        )),
        "reference_mean":a["mean"],
        "candidate_mean":b["mean"],
        "reference_covariance":a["covariance"],
        "candidate_covariance":b["covariance"],
    }


def gaussian_wavefunction_on_grid(coefficients, basis, points):
    return reconstruct_spinor_complete_wavefunction(
        coefficients,basis,points
    )


def compare_wavefunctions(reference, candidate, points, area):
    """Full-wavefunction comparison after removing only a global phase."""
    reference=np.asarray(reference,dtype=complex)
    candidate=np.asarray(candidate,dtype=complex)

    aligned,phase=phase_align_wavefunction(
        reference,candidate,area
    )
    fidelity=phase_aligned_fidelity(
        reference,candidate,area
    )
    l2=phase_aligned_l2_error(
        reference,candidate,area
    )
    density_l2=nuclear_density_l2_error(
        reference,candidate,area
    )
    density_tv=nuclear_density_total_variation(
        reference,candidate,area
    )
    moments=moment_errors(
        reference,candidate,points,area
    )

    return {
        "fidelity":fidelity,
        "phase_aligned_l2_error":l2,
        "alignment_factor":phase,
        "nuclear_density_l2_error":density_l2,
        "nuclear_density_total_variation":density_tv,
        "mean_error_l2":moments["mean_l2"],
        "covariance_error_frobenius":
            moments["covariance_frobenius"],
        "reference_mean":moments["reference_mean"],
        "candidate_mean":moments["candidate_mean"],
        "reference_covariance":
            moments["reference_covariance"],
        "candidate_covariance":
            moments["candidate_covariance"],
        "reference_grid_norm":
            grid_wavefunction_norm(reference,area),
        "candidate_grid_norm":
            grid_wavefunction_norm(candidate,area),
        "aligned_candidate_grid_norm":
            grid_wavefunction_norm(aligned,area),
    }
