from pathlib import Path
import json

root=Path(__file__).resolve().parents[1]
data=json.loads(
    (root/"results"/"v020_sparse_molecular_campaign.json").read_text()
)
idx=data["canonical"]["sparse_provider"]["spatial_index"]

print("v0.20 indexed electronic cache")
print("------------------------------")
for key,value in idx.items():
    print(f"{key}: {value}")
