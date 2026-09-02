from pathlib import Path
import json

root=Path(__file__).resolve().parents[1]
data=json.loads(
    (root/"results"/"v020_sparse_molecular_campaign.json").read_text()
)
c=data["canonical"]

print("v0.20 canonical sparse molecular propagation")
print("--------------------------------------------")
print("possible pairs:",c["total_offdiagonal_pairs"])
print("active edges:",c["final_active_edges"])
print("average sparsity:",c["average_sparsity_fraction"])
print("dense-metric coefficient error:",
      c["metric_coefficient_error"])
print("norm drift:",c["maximum_norm_drift"])
print("final sentinel:",c["sentinels"]["final"])
