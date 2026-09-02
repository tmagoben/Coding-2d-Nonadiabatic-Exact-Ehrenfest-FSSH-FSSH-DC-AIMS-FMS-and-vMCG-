from pathlib import Path
import json

root=Path(__file__).resolve().parents[1]
data=json.loads(
    (root/"results"/"v018_convergence_complete_campaign.json").read_text()
)
dt=data["dt_axis"]

print("v0.18 timestep self convergence")
print("-------------------------------")
for row in dt["rows"]:
    print(
        f"dt={row['dt']:.4f}  "
        f"exact-projected L2="
        f"{row['wavefunction_projected']['phase_aligned_l2_error']:.9f}  "
        f"resolved={row['resolved_control_steps']}"
    )

print("\nSuccessive Gaussian solution differences")
for row in dt["successive_solution_differences"]:
    print(row)

print("\nObserved order:",dt["observed_self_order"])
