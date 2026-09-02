from pathlib import Path
import json

root=Path(__file__).resolve().parents[1]
data=json.loads(
    (root/"results"/"v013_residual_driven_campaign.json").read_text()
)

v12=data["v12_context"]
v13=data["reference"]

print("v0.12 -> v0.13 comparison")
print("-------------------------")
for label,key in [
    ("projection fidelity","projection_fidelity"),
    ("relative residual","relative_residual"),
    ("initial density error","initial_density_error"),
    ("projected dynamics error","projected_dynamics_density_error"),
    ("target density error","target_density_error"),
    ("target population error","target_population_error"),
    ("coherence phase error","coherence_phase_error"),
    ("maximum condition number","max_condition_number"),
]:
    print(
        f"{label:28s} "
        f"v0.12={v12[key]:.9g}   "
        f"v0.13={v13[key]:.9g}"
    )

print("\nv0.13 acceptance:")
print(data["acceptance"])
