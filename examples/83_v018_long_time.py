from pathlib import Path
import json

root=Path(__file__).resolve().parents[1]
data=json.loads(
    (root/"results"/"v018_convergence_complete_campaign.json").read_text()
)

print("v0.18 long-time full-wavefunction trajectory")
print("--------------------------------------------")
for row in data["long_time"]["trajectory"]:
    print(
        f"t={row['time']:.1f} "
        f"fidelity={row['fidelity']:.8f} "
        f"L2={row['phase_aligned_l2_error']:.8f} "
        f"densityL2={row['nuclear_density_l2_error']:.8f}"
    )
