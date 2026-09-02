from pathlib import Path
import json

root=Path(__file__).resolve().parents[1]
data=json.loads(
    (root/"results"/"v014_time_adaptive_defect_campaign.json").read_text()
)

r=data["reference"]
events=[
    e for e in data["adaptive"]["events"]
    if e["kind"]=="defect_enrichment"
]

print("v0.14 time-adaptive TDSE-defect control")
print("---------------------------------------")
print("initial basis:",r["initial_basis_size"])
print("final basis:",r["final_basis_size"])
print("average basis:",r["average_basis_size"])

for e in events:
    print("\nEnrichment")
    print("  step:",e["step"])
    print("  time:",e["time"])
    print("  candidate:",e["candidate_label"])
    print("  candidate count:",e["candidate_count"])
    print("  defect before:",e["relative_defect_before"])
    print("  defect after:",e["relative_defect_after"])
    print("  predicted capture:",e["capture_fraction_predicted"])
    print("  zero coefficient insertion:",e["zero_coefficient_insertion"])

print("\nFinal target-density error:",r["target_density_error"])
print("Projected-state dynamics error:",r["projected_dynamics_density_error"])
