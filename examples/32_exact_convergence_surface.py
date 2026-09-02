import numpy as np

from gaussian_dynamics import (
    CIPassageConfig,
    run_exact_grid_timestep_surface,
    select_finest_exact_reference,
)

rows=run_exact_grid_timestep_surface(
    CIPassageConfig(),
    grid_values=(32,48,64),
    dt_values=(0.010,0.005,0.0025),
)

ref=select_finest_exact_reference(rows)
pref=np.asarray(ref["populations"])

print("Exact grid x timestep surface")
print("-----------------------------")
for row in sorted(rows,key=lambda r:(r["grid_n"],r["dt"])):
    err=np.linalg.norm(np.asarray(row["populations"])-pref)
    print(
        f"N={row['grid_n']:3d}  dt={row['dt']:8.4g}  "
        f"P={np.asarray(row['populations'])}  "
        f"|P-P_finest|={err:.6e}  norm={row['norm']:.15f}"
    )

print("\nFinest candidate:",ref)
print(
    "The finest row is only a candidate reference; inspect its neighboring "
    "refinement differences before calling it converged."
)
