from pathlib import Path
import json

root=Path(__file__).resolve().parents[1]
data=json.loads(
    (root/"results"/"v018_convergence_complete_campaign.json").read_text()
)
c=data["canonical"]["complexity"]

print("v0.18 batched residual-candidate ranking")
print("---------------------------------------")
for key in [
    "candidate_searches",
    "candidates_scored",
    "candidate_batches",
    "candidate_max_dense_grid_elements",
    "candidate_peak_grid_elements",
    "candidate_peak_memory_reduction_fraction",
]:
    print(f"{key:42s}: {c[key]}")
