from pathlib import Path
import json

root=Path(__file__).resolve().parents[1]
data=json.loads(
    (root/"results"/"v017_sparse_error_control_campaign.json").read_text()
)

print("v0.17 bounded-locality S/H/T scoring benchmark")
print("----------------------------------------------")
for row in data["scaling"]:
    print(
        f"N={row['n_basis']:3d} "
        f"active={row['active_edges']:4d} "
        f"exact_checks={row['exact_pair_checks']:4d} "
        f"pair_reduction={100*row['pair_reduction_fraction']:.2f}% "
        f"assembly_speedup={row['assembly_speedup_vs_dense']:.2f}x"
    )

print("\nFitted exponents")
print(data["scaling_fit"])
