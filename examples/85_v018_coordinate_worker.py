from gaussian_dynamics import (
    ConvergenceCoordinatesV18,
    run_coordinate_worker_v18,
)

coordinate=ConvergenceCoordinatesV18(
    dt=0.005,
    max_basis=13,
    local_score_budget=0.01,
    enrich_threshold=0.015,
)

result=run_coordinate_worker_v18(
    coordinate,
    trajectory=False,
)

row=result["result"]
print("Fresh-process v0.18 convergence coordinate")
print("------------------------------------------")
print("coordinates:",row["coordinates"])
print("basis size:",row["basis_size"])
print("projected fidelity:",
      row["projected_wavefunction_fidelity"])
print("projected L2:",
      row["projected_wavefunction_l2_error"])
print("wall seconds:",row["wall_seconds"])
