from pathlib import Path
import json

root=Path(__file__).resolve().parents[1]
data=json.loads(
    (root/"results"/"v012_representation_consistent_campaign.json").read_text()
)

print("v0.12 initial Gaussian-bank projection ladder")
print("----------------------------------------------")
for row in data["projection_ladder"]:
    print(
        f"{row['bank']:>5s}  "
        f"N={row['n_gaussians']:2d}  "
        f"Fproj={row['initial_projection_fidelity']:.8f}  "
        f"init_rho={row['initial_density_error']:.6e}  "
        f"dyn_rho={row['projected_dynamics_density_error']:.6e}  "
        f"target_rho={row['target_density_error']:.6e}  "
        f"pop={row['target_population_error']:.6e}  "
        f"cond={row['max_condition_number']:.3e}"
    )

print("\nAcceptance:")
print(data["acceptance"])
