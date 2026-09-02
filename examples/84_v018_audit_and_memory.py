from pathlib import Path
import json

root=Path(__file__).resolve().parents[1]
data=json.loads(
    (root/"results"/"v018_convergence_complete_campaign.json").read_text()
)
c=data["canonical"]["complexity"]

print("v0.18 sampled-audit and candidate-memory diagnostics")
print("---------------------------------------------------")
print("dense sentinels:",c["sentinel_dense_audits"])
print("sampled audits:",c["sampled_audits"])
print("sampled audit failures:",c["sampled_audit_failures"])
print("v0.18 dense sentinel pair factorizations:",
      c["sentinel_pair_factorizations"])
print("v0.17 dense audit pair factorizations:",
      data["v17_context"]["complexity"]["audit_pair_factorizations"])
print("dense audit work reduction:",
      data["acceptance"]["dense_audit_pair_reduction_vs_v17"])
print("unbatched candidate-grid elements:",
      c["candidate_max_dense_grid_elements"])
print("batched peak elements:",
      c["candidate_peak_grid_elements"])
print("peak memory reduction:",
      c["candidate_peak_memory_reduction_fraction"])

print("\nSampled-audit scaling")
for row in data["sampled_audit_scaling"]:
    print(row)
