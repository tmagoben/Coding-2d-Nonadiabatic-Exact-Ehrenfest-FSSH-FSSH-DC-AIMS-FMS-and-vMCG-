import numpy as np

from gaussian_dynamics import (
    CIPassageConfig,
    AnalyticCI2DFrameProvider,
    DynamicGraphTBF,
    run_basis_complete_graph_aims,
    reduced_electronic_density_analytic_ci_diabatic,
    density_matrix_populations,
    density_matrix_purity,
    basis_completeness_report,
)

config=CIPassageConfig()
provider=AnalyticCI2DFrameProvider(nuclear_mass_au=config.mass)

seed=DynamicGraphTBF(
    uid=0,
    state=config.state,
    q=config.q_array(),
    p=config.p_array(),
    A=config.A_matrix(),
    node=("seed",0),
)

out=run_basis_complete_graph_aims(
    [seed],
    [1.0+0j],
    provider=provider,
    dt=0.005,
    steps=120,
    spa_order=1,
    spawn_action_threshold=1e-4,
    overlap_block=0.9999,
    child_overlap_block=0.995,
    max_basis=10,
    max_generation=5,
    children_per_event=2,
    minimum_spawn_separation_steps=4,
    position_shifts=(0.0,0.05,-0.05),
    width_scales=(0.65,1.0,1.55),
    store_every=20,
)

rho=reduced_electronic_density_analytic_ci_diabatic(
    out["final_coefficients"],
    out["final_basis"],
    normalize=True,
)

report=basis_completeness_report(out)

print("v0.11 strong-CI basis-completeness run")
print("--------------------------------------")
print("final basis size:",len(out["final_basis"]))
print("lineage depth:",report["lineage_depth"])
print("generation histogram:",report["generation_histogram"])
print("width diversity ratio:",report["width_diversity_ratio"])
print("canonical participation ratio:",
      report["canonical_participation_ratio"])
print("diabatic reduced populations:",
      density_matrix_populations(rho))
print("reduced-state purity:",
      density_matrix_purity(rho))
print("final generalized norm:",
      out["records"][-1]["norm"])
print("maximum recorded condition number:",
      max(r["condition_number"] for r in out["records"]))
