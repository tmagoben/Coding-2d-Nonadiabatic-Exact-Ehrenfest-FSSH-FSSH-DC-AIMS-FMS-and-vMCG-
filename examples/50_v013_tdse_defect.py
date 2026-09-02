from pathlib import Path
import json

root=Path(__file__).resolve().parents[1]
data=json.loads(
    (root/"results"/"v013_residual_driven_campaign.json").read_text()
)
d=data["defect_enrichment"]

print("v0.13 instantaneous TDSE-defect enrichment")
print("-------------------------------------------")
print("selected candidate:",d["selected_label"])
print("defect norm before:",d["defect_norm_before"])
print("defect norm after:",d["defect_norm_after"])
print("relative defect before:",d["relative_defect_before"])
print("relative defect after:",d["relative_defect_after"])
print("capture fraction:",d["capture_fraction"])
print("predicted squared reduction:",d["predicted_squared_reduction"])
print("actual squared reduction:",d["actual_squared_reduction"])
print("expanded condition number:",d["expanded_condition_number"])
print("zero-coefficient insertion:",d["zero_coefficient_insertion"])

print(
    "\nThe new Gaussian pair changes the available Galerkin tangent space while "
    "leaving the instantaneous wavefunction unchanged."
)
