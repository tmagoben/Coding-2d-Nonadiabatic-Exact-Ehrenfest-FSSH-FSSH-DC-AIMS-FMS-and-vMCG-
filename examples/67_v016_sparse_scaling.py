from pathlib import Path
import json

root=Path(__file__).resolve().parents[1]
data=json.loads(
    (root/"results"/"v016_sparse_locality_campaign.json").read_text()
)

print("v0.16 bounded-locality scaling benchmark")
print("----------------------------------------")
for row in data["sparse_scaling"]:
    print(
        f"N={row['n_basis']:3d}  "
        f"edges={row['active_edges']:4d}  "
        f"edge_fraction={row['edge_fraction']:.5f}  "
        f"pair_reduction={100*row['pair_reduction_fraction']:.2f}%  "
        f"assembly_speedup={row['assembly_speedup_vs_dense']:.2f}x"
    )

print("\nFitted exponents")
for key,value in data["sparse_scaling_fit"].items():
    print(f"{key}: {value}")
