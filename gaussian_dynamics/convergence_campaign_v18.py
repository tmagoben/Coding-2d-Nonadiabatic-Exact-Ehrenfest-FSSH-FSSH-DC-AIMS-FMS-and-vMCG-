from dataclasses import dataclass, asdict
import math
import numpy as np

from .sparse_pair_matrices_v16 import sparse_reduced_density
from .electronic_observables import (
    density_matrix_populations,
    density_matrix_purity,
)
from .coherence_metrics import (
    coherence_phase_error,
    density_trace_distance,
)
from .wavefunction_metrics_v18 import (
    gaussian_wavefunction_on_grid,
    compare_wavefunctions,
)
from .residual_basis_v13 import normalized_grid_density


@dataclass(frozen=True)
class ConvergenceCoordinatesV18:
    dt: float
    max_basis: int
    local_score_budget: float
    enrich_threshold: float

    def as_dict(self):
        return asdict(self)


def _reduced_density_metrics(rho,reference):
    rho=np.asarray(rho,dtype=complex)
    ref=np.asarray(reference,dtype=complex)
    return {
        "density_frobenius_error":float(
            np.linalg.norm(rho-ref,ord="fro")
        ),
        "trace_distance":
            density_trace_distance(rho,ref),
        "population_l2_error":float(
            np.linalg.norm(
                density_matrix_populations(rho)
                -density_matrix_populations(ref)
            )
        ),
        "purity":density_matrix_purity(rho),
        "reference_purity":
            density_matrix_purity(ref),
        "purity_error":float(abs(
            density_matrix_purity(rho)
            -density_matrix_purity(ref)
        )),
        "coherence_phase_error":
            coherence_phase_error(rho,ref),
    }


def evaluate_convergence_run_v18(
    run_output,
    comparison_grid,
    exact_projected_final,
    exact_target_final,
):
    """Evaluate reduced-state and full-wavefunction errors at one final time."""
    rho=sparse_reduced_density(
        run_output["final_coefficients"],
        run_output["final_sparse_matrices"].Snuc,
        normalize=True,
    )
    psi=gaussian_wavefunction_on_grid(
        run_output["final_coefficients"],
        run_output["final_basis"],
        comparison_grid.points,
    )
    area=float(comparison_grid.dx*comparison_grid.dx)

    rho_projected=normalized_grid_density(
        exact_projected_final,
        comparison_grid.dx,
    )
    rho_target=normalized_grid_density(
        exact_target_final,
        comparison_grid.dx,
    )

    return {
        "basis_size":int(len(run_output["final_basis"])),
        "average_basis_size":
            float(run_output["average_basis_size"]),
        "reduced_density_projected":
            _reduced_density_metrics(
                rho,rho_projected
            ),
        "reduced_density_target":
            _reduced_density_metrics(
                rho,rho_target
            ),
        "wavefunction_projected":
            compare_wavefunctions(
                exact_projected_final,
                psi,
                comparison_grid.points,
                area,
            ),
        "wavefunction_target":
            compare_wavefunctions(
                exact_target_final,
                psi,
                comparison_grid.points,
                area,
            ),
        "maximum_norm_drift":float(max(
            abs(row["norm"]-1.0)
            for row in run_output["records"]
        )),
        "maximum_condition_number":float(max(
            row["condition_number"]
            for row in run_output["records"]
        )),
        "final_density_matrix":rho,
        "complexity":
            run_output["complexity"],
    }


def compare_snapshot_trajectory_v18(
    snapshots,
    exact_run,
    comparison_grid,
):
    """Full-wavefunction metrics at matching stored times."""
    exact_times=np.asarray(exact_run["time"],dtype=float)
    exact_states=np.asarray(exact_run["psi"],dtype=complex)
    area=float(comparison_grid.dx*comparison_grid.dx)

    rows=[]
    for snap in snapshots:
        t=float(snap["time"])
        k=int(np.argmin(np.abs(exact_times-t)))
        if abs(exact_times[k]-t)>1e-9:
            raise ValueError(
                f"no exact reference stored at t={t}"
            )

        psi_g=gaussian_wavefunction_on_grid(
            snap["coefficients"],
            snap["basis"],
            comparison_grid.points,
        )
        metrics=compare_wavefunctions(
            exact_states[k],
            psi_g,
            comparison_grid.points,
            area,
        )
        rows.append({
            "step":int(snap["step"]),
            "time":t,
            **metrics,
        })
    return rows


def observed_order_from_dt(rows,error_key):
    """Pairwise observed order for rows sorted coarse -> fine in dt."""
    rows=sorted(rows,key=lambda x:x["coordinates"]["dt"],reverse=True)
    out=[]
    for a,b in zip(rows[:-1],rows[1:]):
        h1=float(a["coordinates"]["dt"])
        h2=float(b["coordinates"]["dt"])
        e1=float(a[error_key])
        e2=float(b[error_key])
        if e1<=0.0 or e2<=0.0 or h1==h2:
            p=None
        else:
            p=float(math.log(e1/e2)/math.log(h1/h2))
        out.append({
            "dt_coarse":h1,
            "dt_fine":h2,
            "error_coarse":e1,
            "error_fine":e2,
            "observed_order":p,
        })
    return out


def axis_sensitivity_summary(rows,metric_key,higher_is_better=False):
    values=np.asarray([
        float(row[metric_key])
        for row in rows
    ],dtype=float)
    if higher_is_better:
        best=int(np.argmax(values))
        worst=int(np.argmin(values))
    else:
        best=int(np.argmin(values))
        worst=int(np.argmax(values))
    return {
        "minimum":float(np.min(values)),
        "maximum":float(np.max(values)),
        "span":float(np.max(values)-np.min(values)),
        "higher_is_better":bool(higher_is_better),
        "best_index":best,
        "worst_index":worst,
    }


def refinement_ladder_summary(rows,error_key):
    if len(rows)<2:
        raise ValueError("refinement ladder requires at least two rows.")
    coarse=float(rows[0][error_key])
    fine=float(rows[-1][error_key])
    return {
        "coarse_error":coarse,
        "fine_error":fine,
        "absolute_improvement":float(coarse-fine),
        "relative_improvement_fraction":float(
            (coarse-fine)/max(abs(coarse),1e-30)
        ),
        "fine_better_than_coarse":bool(fine<coarse),
        "strictly_monotone":bool(all(
            float(b[error_key])<float(a[error_key])
            for a,b in zip(rows[:-1],rows[1:])
        )),
    }


def successive_self_convergence_order(
    coarse_medium_error,
    medium_fine_error,
    refinement_ratio=2.0,
):
    """Observed order from successive solution differences.

    For an asymptotic method with error C h^p,

        ||u_h-u_{h/r}|| / ||u_{h/r}-u_{h/r^2}|| -> r^p.

    This avoids treating the finest numerical solution as an exact reference.
    """
    e1=float(coarse_medium_error)
    e2=float(medium_fine_error)
    r=float(refinement_ratio)
    if e1<=0.0 or e2<=0.0 or r<=1.0:
        return None
    return float(np.log(e1/e2)/np.log(r))
