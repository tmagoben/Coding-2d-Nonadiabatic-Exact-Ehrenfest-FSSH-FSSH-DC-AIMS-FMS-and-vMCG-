from pathlib import Path
import json

root=Path(__file__).resolve().parents[1]
data=json.loads(
    (root/"results"/"v019_molecular_direct_dynamics_campaign.json").read_text()
)
p=data["provider_scan"]

print("v0.19 Cartesian molecular -> generalized projection")
print("--------------------------------------------------")
print("scan points:",p["points"])
print("tracked max energy error:",
      p["maximum_tracked_energy_error"])
print("tracked max gradient error:",
      p["maximum_tracked_gradient_error"])
print("tracked max NAC error:",
      p["maximum_tracked_nac_error"])
print("raw scrambled max energy error:",
      p["raw_scrambled_max_energy_error"])
