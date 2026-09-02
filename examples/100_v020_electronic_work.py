from pathlib import Path
import json

root=Path(__file__).resolve().parents[1]
data=json.loads(
    (root/"results"/"v020_sparse_molecular_campaign.json").read_text()
)
c=data["canonical"]

print("v0.20 molecular electronic-work comparison")
print("------------------------------------------")
print("sparse cache misses:",
      c["sparse_provider"]["cache_misses"])
print("dense cache misses:",
      c["dense_provider"]["cache_misses"])
print("reduction:",
      c["backend_miss_reduction_fraction"])
print("diagnostic wall speedup:",
      c["wall_speedup"])
