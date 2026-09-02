from pathlib import Path
import json

root=Path(__file__).resolve().parents[1]
data=json.loads(
    (root/"results"/"v020_sparse_molecular_campaign.json").read_text()
)

print("v0.20 bounded-locality molecular scaling")
print("----------------------------------------")
for row in data["scaling"]:
    print(
        f"N={row['n_basis']:3d}  "
        f"active={row['active_edges']:4d}  "
        f"checks={row['exact_pair_checks']:4d}  "
        f"dense={row['dense_offdiagonal_pairs']:6d}  "
        f"reduction={100*row['pair_check_reduction_fraction']:.2f}%  "
        f"ES={row['provider_cache_misses']:4d}"
    )

print("\nFitted exponents")
for key,value in data["scaling_fit"].items():
    print(f"{key}: {value}")
