from pathlib import Path
import json

root=Path(__file__).resolve().parents[1]
data=json.loads(
    (root/"results"/"v015_cost_aware_cache_campaign.json").read_text()
)

event=[
    e for e in data["adaptive"]["events"]
    if e["kind"]=="cost_aware_defect_enrichment"
][0]

print("v0.15 cost-aware TDSE-defect event")
print("----------------------------------")
print("step:",event["step"])
print("time:",event["time"])
print("candidate:",event["candidate_label"])
print("defect before:",event["relative_defect_before"])
print("defect after:",event["relative_defect_after"])
print("capture fraction:",event["capture_fraction_predicted"])
print("normalized incremental cost:",event["normalized_incremental_cost"])
print("cost-aware utility:",event["cost_aware_utility"])
print("estimated horizon seconds:",event["estimated_incremental_seconds"])
print("expanded condition:",event["expanded_condition_number"])
print(
    "new pair factorizations during accepted matrix expansion:",
    event["new_pair_factorizations_during_expansion"],
)
