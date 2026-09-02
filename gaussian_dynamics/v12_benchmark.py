from dataclasses import dataclass, asdict
from pathlib import Path
import json
import numpy as np

from .benchmark_campaign import CIPassageConfig
from .dynamic_gauge_graph import AnalyticCI2DFrameProvider
from .dynamic_graph_aims import DynamicGraphTBF
from .born_huang_grid_v12 import build_born_huang_grid_2d
from .exact_benchmark import localized_adiabatic_packet_2d
from .exact2d import run_exact_2d
from .ci2d import diabatic_potential_2d
from .electronic_observables import (
    exact_reduced_electronic_density_diabatic,
    density_matrix_populations,
    density_matrix_purity,
    density_matrix_linear_entropy,
    density_matrix_von_neumann_entropy,
)
from .coherence_metrics import (
    coherence_magnitude,
    coherence_magnitude_error,
    coherence_phase_error,
    density_trace_distance,
    bloch_vector,
    bloch_vector_error,
)
from .initial_projection_v12 import (
    project_grid_wavefunction_to_spinor_complete_basis,
)
from .spinor_complete_dynamics_v12 import (
    run_spinor_complete_lvc_gaussians,
)
from .spinor_complete_lvc_v12 import (
    spinor_complete_reduced_density,
)


@dataclass(frozen=True)
class V12AcceptanceThresholds:
    max_initial_density_error: float = 0.05
    max_projected_dynamics_density_error: float = 1e-3
    max_target_density_error: float = 0.05
    max_target_population_error: float = 0.05
    max_coherence_phase_error: float = 0.01
    max_norm_drift: float = 1e-4
    max_condition_number: float = 1e5


def _normalized_density(psi,dx):
    rho=exact_reduced_electronic_density_diabatic(
        psi,dx,dx
    )
    return rho/np.trace(rho)


def _bank(config,kind):
    q0=config.q_array()
    p0=config.p_array()
    A0=config.A_matrix()

    if kind=="one":
        shifts=[(0.0,0.0)]
        scale=1.0
    elif kind=="five":
        shifts=[
            (0.0,0.0),
            (0.35,0.0),(-0.35,0.0),
            (0.0,0.35),(0.0,-0.35),
        ]
        scale=2.0
    elif kind=="nine":
        shifts=[
            (i*0.45,j*0.45)
            for i in (-1,0,1)
            for j in (-1,0,1)
        ]
        scale=3.0
    else:
        raise ValueError("kind must be one, five, or nine.")

    basis=[]
    for k,shift in enumerate(shifts):
        basis.append(
            DynamicGraphTBF(
                uid=k,
                state=config.state,
                q=q0+np.asarray(shift,float),
                p=p0.copy(),
                A=scale*A0,
                node=("v12_initial_bank",kind,k),
            )
        )
    return basis


def _case_metrics(
    config,
    provider,
    grid,
    psi_target,
    rho_target_initial,
    rho_target_final,
    V,
    bank_kind,
):
    basis=_bank(config,bank_kind)

    projection=project_grid_wavefunction_to_spinor_complete_basis(
        psi_target,
        grid.points,
        grid.dx,
        basis,
        provider,
    )

    rho_projected_initial=_normalized_density(
        projection.projected_wavefunction,
        grid.dx,
    )

    exact_projected=run_exact_2d(
        projection.projected_wavefunction,
        grid.dx,
        grid.dx,
        V,
        mass=config.mass,
        dt=0.0025,
        steps=int(round(config.final_time/0.0025)),
        store_every=int(round(config.final_time/0.0025)),
    )
    rho_exact_projected_final=_normalized_density(
        exact_projected["psi"][-1],
        grid.dx,
    )

    gaussian=run_spinor_complete_lvc_gaussians(
        basis,
        C0=projection.coefficients,
        provider=provider,
        dt=0.005,
        steps=int(round(config.final_time/0.005)),
        integrator="cayley",
        spawn_action_threshold=1e9,
        max_basis=len(basis),
        condition_limit=1e12,
        store_every=max(1,int(round(config.final_time/0.005))//6),
    )

    rho_gaussian=spinor_complete_reduced_density(
        gaussian["final_coefficients"],
        gaussian["final_basis"],
        normalize=True,
    )

    target_pop=density_matrix_populations(rho_target_final)
    gauss_pop=density_matrix_populations(rho_gaussian)

    return {
        "bank":bank_kind,
        "n_gaussians":len(basis),
        "initial_projection_fidelity":float(projection.fidelity),
        "initial_projection_relative_residual":float(projection.relative_residual),
        "initial_projection_condition_number":float(projection.condition_number),
        "initial_density_error":float(np.linalg.norm(
            rho_projected_initial-rho_target_initial,
            ord="fro",
        )),
        "projected_exact_to_target_density_error":float(np.linalg.norm(
            rho_exact_projected_final-rho_target_final,
            ord="fro",
        )),
        "projected_dynamics_density_error":float(np.linalg.norm(
            rho_gaussian-rho_exact_projected_final,
            ord="fro",
        )),
        "target_density_error":float(np.linalg.norm(
            rho_gaussian-rho_target_final,
            ord="fro",
        )),
        "target_trace_distance":density_trace_distance(
            rho_gaussian,rho_target_final
        ),
        "target_population_error":float(np.linalg.norm(
            gauss_pop-target_pop
        )),
        "target_populations":target_pop,
        "gaussian_populations":gauss_pop,
        "target_purity":density_matrix_purity(rho_target_final),
        "gaussian_purity":density_matrix_purity(rho_gaussian),
        "purity_error":float(abs(
            density_matrix_purity(rho_gaussian)
            -density_matrix_purity(rho_target_final)
        )),
        "target_coherence":complex(rho_target_final[0,1]),
        "gaussian_coherence":complex(rho_gaussian[0,1]),
        "coherence_magnitude_error":coherence_magnitude_error(
            rho_gaussian,rho_target_final
        ),
        "coherence_phase_error":coherence_phase_error(
            rho_gaussian,rho_target_final
        ),
        "target_bloch_vector":bloch_vector(rho_target_final),
        "gaussian_bloch_vector":bloch_vector(rho_gaussian),
        "bloch_vector_error":bloch_vector_error(
            rho_gaussian,rho_target_final
        ),
        "max_norm_drift":float(max(
            abs(r["norm"]-1.0)
            for r in gaussian["records"]
        )),
        "max_condition_number":float(max(
            r["condition_number_nuclear"]
            for r in gaussian["records"]
        )),
    }


def evaluate_v12_acceptance(reference_case,thresholds=None):
    t=thresholds or V12AcceptanceThresholds()
    phase=reference_case["coherence_phase_error"]
    checks={
        "initial_density_representation":
            reference_case["initial_density_error"] <= t.max_initial_density_error,
        "projected_dynamics":
            reference_case["projected_dynamics_density_error"]
            <= t.max_projected_dynamics_density_error,
        "target_full_density":
            reference_case["target_density_error"] <= t.max_target_density_error,
        "target_population":
            reference_case["target_population_error"] <= t.max_target_population_error,
        "coherence_phase":
            phase is not None and phase <= t.max_coherence_phase_error,
        "norm":
            reference_case["max_norm_drift"] <= t.max_norm_drift,
        "conditioning":
            reference_case["max_condition_number"] <= t.max_condition_number,
    }
    return {
        "passed":bool(all(checks.values())),
        "checks":checks,
        "thresholds":asdict(t),
    }


def _decode_complex_array(value):
    arr=np.asarray(value)
    if arr.ndim>=1 and arr.shape[-1]==2:
        return arr[...,0].astype(float)+1j*arr[...,1].astype(float)
    return np.asarray(value,dtype=complex)


def load_v11_release_context(repository_root):
    path=Path(repository_root)/"results"/"v011_basis_completeness_campaign.json"
    if not path.exists():
        return None
    data=json.loads(path.read_text(encoding="utf-8"))

    rho_v11=_decode_complex_array(data["v11_reference_rho"])
    rho_exact=_decode_complex_array(data["exact"]["rho"])

    return {
        "v11_population_error":data["v11_reference"]["population_l2_error"],
        "v11_center_frozen_density_error":data["v11_reference"]["density_frobenius_error"],
        "v11_trace_distance":density_trace_distance(rho_v11,rho_exact),
        "v11_purity":data["v11_reference"]["purity"],
        "v11_purity_error":float(abs(
            density_matrix_purity(rho_v11)-density_matrix_purity(rho_exact)
        )),
        "v11_coherence_magnitude":coherence_magnitude(rho_v11),
        "exact_coherence_magnitude":coherence_magnitude(rho_exact),
        "v11_coherence_phase_error":coherence_phase_error(rho_v11,rho_exact),
        "v11_bloch_vector_error":bloch_vector_error(rho_v11,rho_exact),
        "v11_acceptance":data["acceptance"],
    }


def run_v012_release_benchmark(
    config=CIPassageConfig(),
    repository_root=None,
):
    provider=AnalyticCI2DFrameProvider(
        nuclear_mass_au=config.mass
    )
    grid=build_born_huang_grid_2d(
        grid_n=64,
        half_width=config.half_width,
        mass=config.mass,
        params=provider.params,
    )

    psi_target=localized_adiabatic_packet_2d(
        grid.points,
        config.q_array(),
        config.p_array(),
        config.A_matrix(),
        state=config.state,
    )

    rho_target_initial=_normalized_density(
        psi_target,grid.dx
    )

    V=diabatic_potential_2d(
        grid.X,grid.Y,provider.params
    )

    exact_target=run_exact_2d(
        psi_target,
        grid.dx,
        grid.dx,
        V,
        mass=config.mass,
        dt=0.0025,
        steps=int(round(config.final_time/0.0025)),
        store_every=int(round(config.final_time/0.0025)),
    )
    rho_target_final=_normalized_density(
        exact_target["psi"][-1],
        grid.dx,
    )

    cases=[
        _case_metrics(
            config,provider,grid,psi_target,
            rho_target_initial,rho_target_final,V,kind
        )
        for kind in ("one","five","nine")
    ]

    reference=cases[-1]
    acceptance=evaluate_v12_acceptance(reference)

    center_spinor=np.asarray(
        provider.evaluate(config.q_array()).frame[:,int(config.state)],
        dtype=complex,
    )
    rho_center=np.outer(center_spinor,np.conj(center_spinor))
    rho_center/=np.trace(rho_center)

    center_frozen_initial={
        "density":rho_center,
        "density_error":float(np.linalg.norm(
            rho_center-rho_target_initial,
            ord="fro",
        )),
        "populations":density_matrix_populations(rho_center),
        "purity":density_matrix_purity(rho_center),
        "coherence_magnitude":coherence_magnitude(rho_center),
    }

    context=None
    if repository_root is not None:
        context=load_v11_release_context(repository_root)

    return {
        "config":{
            "q0":list(config.q0),
            "p0":list(config.p0),
            "A_diag":list(config.A_diag),
            "state":config.state,
            "mass":config.mass,
            "final_time":config.final_time,
            "half_width":config.half_width,
        },
        "exact_target":{
            "initial_density":rho_target_initial,
            "initial_populations":density_matrix_populations(rho_target_initial),
            "initial_purity":density_matrix_purity(rho_target_initial),
            "initial_coherence_magnitude":coherence_magnitude(rho_target_initial),
            "final_density":rho_target_final,
            "final_populations":density_matrix_populations(rho_target_final),
            "final_purity":density_matrix_purity(rho_target_final),
            "final_linear_entropy":density_matrix_linear_entropy(rho_target_final),
            "final_von_neumann_entropy":
                density_matrix_von_neumann_entropy(rho_target_final),
            "final_coherence_magnitude":coherence_magnitude(rho_target_final),
        },
        "center_frozen_initial":center_frozen_initial,
        "projection_ladder":cases,
        "reference_case":reference,
        "acceptance":acceptance,
        "v11_context":context,
    }
