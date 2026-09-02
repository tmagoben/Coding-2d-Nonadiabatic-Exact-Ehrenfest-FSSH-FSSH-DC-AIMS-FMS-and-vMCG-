from pathlib import Path
import json

root=Path(__file__).resolve().parents[1]
data=json.loads(
    (root/"results"/"v012_representation_consistent_campaign.json").read_text()
)
v11=data["v11_context"]
v12=data["reference_case"]

print("v0.11 -> v0.12 coherence comparison")
print("------------------------------------")
print("v0.11 full-density error:          ",
      v11["v11_center_frozen_density_error"])
print("v0.12 full-density error:          ",
      v12["target_density_error"])
print()
print("v0.11 coherence phase error / rad: ",
      v11["v11_coherence_phase_error"])
print("v0.12 coherence phase error / rad: ",
      v12["coherence_phase_error"])
print()
print("v0.11 population error:            ",
      v11["v11_population_error"])
print("v0.12 population error:            ",
      v12["target_population_error"])

print(
    "\nv0.12 is evaluated with an explicit initial-state representation audit, "
    "so the density/coherence comparison is no longer interpreted as propagation "
    "error alone."
)
