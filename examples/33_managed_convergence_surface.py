from gaussian_dynamics import (
    CIPassageConfig,
    run_managed_parameter_surface,
)

rows=run_managed_parameter_surface(
    CIPassageConfig(final_time=0.18),
    dts=(0.010,0.005),
    spa_orders=(0,1),
    spawn_action_thresholds=(2e-4,),
    max_basis_values=(2,4),
    overlap_blocks=(0.90,0.9999),
)

print("Managed graph-Gaussian convergence surface")
print("------------------------------------------")
for row in rows:
    print(
        f"dt={row['dt']:7.4f}  SPA={row['spa_order']}  "
        f"Nmax={row['max_basis']}  block={row['overlap_block']:.4f}  "
        f"basis_used={row['max_basis_size']}  spawns={row['spawn_count']}  "
        f"cond_max={row['max_condition_number']:.4e}  "
        f"norm_err={row['max_norm_error']:.3e}"
    )

print(
    "\nThe stored state-label populations are an internal convergence proxy. "
    "Use example 31 for a rigorous exact/Gaussian population comparison."
)
