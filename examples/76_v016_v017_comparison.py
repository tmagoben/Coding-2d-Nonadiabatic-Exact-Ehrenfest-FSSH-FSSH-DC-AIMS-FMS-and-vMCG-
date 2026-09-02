from pathlib import Path
import json

root=Path(__file__).resolve().parents[1]
data=json.loads(
    (root/"results"/"v017_sparse_error_control_campaign.json").read_text()
)

old=data["v16_context"]["reference"]
new=data["reference"]

print("v0.16 -> v0.17 physical comparison")
print("----------------------------------")
for label,key in [
    ("projected dynamics error","projected_dynamics_density_error"),
    ("target density error","target_density_error"),
    ("population error","target_population_error"),
    ("coherence phase error","coherence_phase_error"),
    ("purity","purity"),
]:
    print(
        f"{label:28s} "
        f"v0.16={old[key]:.12g} "
        f"v0.17={new[key]:.12g}"
    )

print(
    "\nfinal density-matrix difference:",
    data["acceptance"]["final_rho_difference_vs_v16"],
)
print(
    "online score relaxations:",
    data["adaptive"]["complexity"]["score_relaxations"],
)
