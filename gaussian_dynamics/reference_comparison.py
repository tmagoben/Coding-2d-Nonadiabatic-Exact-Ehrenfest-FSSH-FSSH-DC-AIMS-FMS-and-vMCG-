import numpy as np

from .benchmark_campaign import (
    CIPassageConfig,
    run_exact_passage,
    run_managed_passage,
)
from .dynamic_gauge_graph import AnalyticCI2DFrameProvider
from .electronic_observables import (
    reduced_electronic_density_graph,
    reduced_electronic_density_analytic_ci_diabatic,
    exact_reduced_electronic_density_diabatic,
    exact_reference_frame_density,
    density_matrix_populations,
    density_matrix_purity,
    density_matrix_linear_entropy,
    density_matrix_von_neumann_entropy,
)


def compare_managed_exact_common_frame(
    config=CIPassageConfig(),
    managed_dt=0.005,
    exact_dt=0.0025,
    exact_grid_n=64,
    spa_order=0,
    spawn_action_threshold=2e-4,
    max_basis=4,
    overlap_block=0.90,
):
    """Compare exact and graph-Gaussian reduced electronic density matrices.

    Both are expressed in the same *fixed* electronic frame located at the final
    center of TBF 0.  This avoids comparing an exact global adiabatic population with
    a basis-label population proxy that is not a projector expectation value.
    """
    exact = run_exact_passage(
        config=config,
        grid_n=exact_grid_n,
        dt=exact_dt,
    )

    managed = run_managed_passage(
        config=config,
        dt=managed_dt,
        spa_order=spa_order,
        spawn_action_threshold=spawn_action_threshold,
        max_basis=max_basis,
        overlap_block=overlap_block,
    )

    reference_tbf = managed["final_basis"][0]
    reference_node = reference_tbf.node

    rho_managed = reduced_electronic_density_graph(
        managed["final_coefficients"],
        managed["final_basis"],
        managed["registry"],
        reference_node,
        normalize=True,
    )

    frame = AnalyticCI2DFrameProvider(
        nuclear_mass_au=config.mass
    ).evaluate(reference_tbf.q).frame

    rho_exact = exact_reference_frame_density(
        exact["psi_final"],
        exact["dx"],
        exact["dx"],
        frame,
        normalize=True,
    )

    return {
        "reference_coordinate": reference_tbf.q.copy(),
        "rho_exact": rho_exact,
        "rho_managed": rho_managed,
        "populations_exact": density_matrix_populations(rho_exact),
        "populations_managed": density_matrix_populations(rho_managed),
        "density_frobenius_error": float(
            np.linalg.norm(rho_managed-rho_exact, ord="fro")
        ),
        "population_l2_error": float(
            np.linalg.norm(
                density_matrix_populations(rho_managed)
                - density_matrix_populations(rho_exact)
            )
        ),
        "purity_exact": density_matrix_purity(rho_exact),
        "purity_managed": density_matrix_purity(rho_managed),
        "linear_entropy_exact": density_matrix_linear_entropy(rho_exact),
        "linear_entropy_managed": density_matrix_linear_entropy(rho_managed),
        "von_neumann_entropy_exact": density_matrix_von_neumann_entropy(rho_exact),
        "von_neumann_entropy_managed": density_matrix_von_neumann_entropy(rho_managed),
        "managed": managed,
        "exact": exact,
    }


def compare_managed_exact_diabatic_density(
    config=CIPassageConfig(),
    managed_dt=0.005,
    exact_dt=0.0025,
    exact_grid_n=64,
    spa_order=0,
    spawn_action_threshold=2e-4,
    max_basis=4,
    overlap_block=0.90,
):
    """Preferred analytic-model comparison in the model's global diabatic basis."""
    exact = run_exact_passage(
        config=config,
        grid_n=exact_grid_n,
        dt=exact_dt,
    )

    managed = run_managed_passage(
        config=config,
        dt=managed_dt,
        spa_order=spa_order,
        spawn_action_threshold=spawn_action_threshold,
        max_basis=max_basis,
        overlap_block=overlap_block,
    )

    rho_exact = exact_reduced_electronic_density_diabatic(
        exact["psi_final"],
        exact["dx"],
        exact["dx"],
    )
    rho_exact = rho_exact/np.trace(rho_exact)

    rho_managed = reduced_electronic_density_analytic_ci_diabatic(
        managed["final_coefficients"],
        managed["final_basis"],
        normalize=True,
    )

    return {
        "rho_exact": rho_exact,
        "rho_managed": rho_managed,
        "populations_exact": density_matrix_populations(rho_exact),
        "populations_managed": density_matrix_populations(rho_managed),
        "density_frobenius_error": float(
            np.linalg.norm(rho_managed-rho_exact, ord="fro")
        ),
        "population_l2_error": float(
            np.linalg.norm(
                density_matrix_populations(rho_managed)
                - density_matrix_populations(rho_exact)
            )
        ),
        "purity_exact": density_matrix_purity(rho_exact),
        "purity_managed": density_matrix_purity(rho_managed),
        "linear_entropy_exact": density_matrix_linear_entropy(rho_exact),
        "linear_entropy_managed": density_matrix_linear_entropy(rho_managed),
        "von_neumann_entropy_exact": density_matrix_von_neumann_entropy(rho_exact),
        "von_neumann_entropy_managed": density_matrix_von_neumann_entropy(rho_managed),
        "managed": managed,
        "exact": exact,
    }
