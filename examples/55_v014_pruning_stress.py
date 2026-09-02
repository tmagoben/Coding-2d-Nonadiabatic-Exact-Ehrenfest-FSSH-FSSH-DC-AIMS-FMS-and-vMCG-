from pathlib import Path
import json

root=Path(__file__).resolve().parents[1]
data=json.loads(
    (root/"results"/"v014_time_adaptive_defect_campaign.json").read_text()
)
p=data["pruning_stress"]

print("v0.14 exact low-loss pruning stress test")
print("----------------------------------------")
print("removed uid:",p["removed_uid"])
print("fractional projection loss:",p["fractional_projection_loss"])
print("condition before:",p["condition_before"])
print("condition after:",p["condition_after"])
print("condition improvement factor:",p["condition_improvement_factor"])

print(
    "\nThe stress Gaussian carries zero electronic amplitude, so its exact "
    "leave-one-out represented-state loss is zero."
)
