from pathlib import Path
import json

root=Path(__file__).resolve().parents[1]
data=json.loads(
    (root/"results"/"v012_representation_consistent_campaign.json").read_text()
)

exact=data["exact_target"]
center=data["center_frozen_initial"]

print("v0.12 initial electronic-representation audit")
print("---------------------------------------------")
print("exact coordinate-dependent initial populations:",
      exact["initial_populations"])
print("exact coordinate-dependent initial purity:",
      exact["initial_purity"])
print()
print("center-frozen initial populations:",
      center["populations"])
print("center-frozen initial purity:",
      center["purity"])
print("center-frozen reduced-density error:",
      center["density_error"])

print(
    "\nThe discrepancy is present at t=0, before any nonadiabatic propagation."
)
