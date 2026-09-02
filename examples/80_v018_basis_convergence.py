from pathlib import Path
import json

root=Path(__file__).resolve().parents[1]
data=json.loads(
    (root/"results"/"v018_convergence_complete_campaign.json").read_text()
)

print("v0.18 basis-completeness ladder")
print("-------------------------------")
for row in data["basis_axis"]:
    print(
        f"Nmax={row['max_basis']:2d}  "
        f"Nfinal={row['final_basis_size']:2d}  "
        f"F={row['wavefunction_projected']['fidelity']:.9f}  "
        f"L2={row['wavefunction_projected']['phase_aligned_l2_error']:.9f}"
    )

print(
    "\nRelative 10 -> 13 improvement:",
    data["acceptance"]["basis_ladder_improvement_fraction"],
)
