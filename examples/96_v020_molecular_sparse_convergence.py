from pathlib import Path
import json

root=Path(__file__).resolve().parents[1]
data=json.loads(
    (root/"results"/"v020_sparse_molecular_campaign.json").read_text()
)
conv=data["molecular_sparse_convergence"]

print("v0.20 molecular score-threshold convergence")
print("-------------------------------------------")
for row in conv["threshold_rows"]:
    print(
        f"eta={row['enter_score']:.3f}  "
        f"E={row['active_edges']:3d}  "
        f"S={row['S_error']:.8g}  "
        f"H={row['H_error']:.8g}  "
        f"T={row['T_error']:.8g}"
    )

print("\nv0.20 local-score-budget convergence")
print("------------------------------------")
for row in conv["budget_rows"]:
    print(
        f"B={row['budget']:.4g}  "
        f"E={row['active_edges']:3d}  "
        f"promoted={row['budget_promoted_edges']:3d}  "
        f"S={row['S_error']:.8g}  "
        f"H={row['H_error']:.8g}  "
        f"T={row['T_error']:.8g}"
    )
