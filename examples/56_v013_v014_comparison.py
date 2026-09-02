from pathlib import Path
import json

root=Path(__file__).resolve().parents[1]
data=json.loads(
    (root/"results"/"v014_time_adaptive_defect_campaign.json").read_text()
)

old=data["v13_context"]
new=data["reference"]

print("v0.13 -> v0.14 comparison")
print("-------------------------")
for label,key in [
    ("projected dynamics error","projected_dynamics_density_error"),
    ("target density error","target_density_error"),
    ("target population error","target_population_error"),
    ("coherence phase error","coherence_phase_error"),
    ("maximum condition","max_condition_number"),
]:
    print(
        f"{label:28s}"
        f" v0.13={old[key]:.9g}"
        f" v0.14={new[key]:.9g}"
    )

print("\nv0.13 basis size:",old["basis_size"])
print("v0.14 initial basis:",new["initial_basis_size"])
print("v0.14 final basis:",new["final_basis_size"])
print("v0.14 average basis:",new["average_basis_size"])

print(
    "\nv0.14 is an adaptive-control upgrade, not a claim that a delayed "
    "11th basis direction must outperform a static 11-Gaussian basis in every "
    "observable."
)
