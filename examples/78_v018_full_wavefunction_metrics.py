from pathlib import Path
import json

root=Path(__file__).resolve().parents[1]
data=json.loads(
    (root/"results"/"v018_convergence_complete_campaign.json").read_text()
)
c=data["canonical"]

print("v0.18 full-wavefunction observables")
print("----------------------------------")
for key in [
    "projected_wavefunction_fidelity",
    "projected_wavefunction_l2_error",
    "projected_nuclear_density_l2_error",
    "projected_nuclear_density_tv",
    "projected_mean_error",
    "projected_covariance_error",
    "projected_reduced_density_error",
]:
    print(f"{key:42s}: {c[key]}")
