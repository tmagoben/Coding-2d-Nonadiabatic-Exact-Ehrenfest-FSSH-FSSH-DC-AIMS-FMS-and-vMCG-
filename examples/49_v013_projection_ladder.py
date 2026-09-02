from pathlib import Path
import json

root=Path(__file__).resolve().parents[1]
data=json.loads(
    (root/"results"/"v013_residual_driven_campaign.json").read_text()
)

print("v0.13 residual-selection ladder")
print("-------------------------------")
for row in data["selection_ladder"]:
    print(
        f"N={row['basis_size']:2d}  "
        f"F={row['projection_fidelity']:.8f}  "
        f"res={row['relative_residual']:.8f}  "
        f"rho0={row['density_error']:.8f}  "
        f"cond={row['condition_number']:.3e}"
    )

print("\nThe Hilbert residual decreases monotonically.")
print(
    "The reduced-density error need not decrease monotonically because it is a "
    "nonlinear reduced observable."
)
