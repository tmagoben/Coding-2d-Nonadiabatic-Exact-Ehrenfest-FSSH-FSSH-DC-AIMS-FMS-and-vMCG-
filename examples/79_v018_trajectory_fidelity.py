from pathlib import Path
import json

root=Path(__file__).resolve().parents[1]
data=json.loads(
    (root/"results"/"v018_convergence_complete_campaign.json").read_text()
)

print("v0.18 projected-reference trajectory")
print("------------------------------------")
for row in data["trajectory_projected_reference"]:
    print(
        f"t={row['time']:.1f}  "
        f"F={row['fidelity']:.12f}  "
        f"L2={row['phase_aligned_l2_error']:.9f}  "
        f"density_L2={row['nuclear_density_l2_error']:.9f}"
    )
