from pathlib import Path
import json

root=Path(__file__).resolve().parents[1]
result=json.loads(
    (root/"results"/"v011_basis_completeness_campaign.json").read_text()
)

print("v0.11 saved branching ablation study")
print("------------------------------------")
print(
    "v0.10 baseline population error:",
    result["v10_baseline"]["population_l2_error"],
)
print(
    "v0.11 reference population error:",
    result["v11_reference"]["population_l2_error"],
)

for name,row in result["ablations"].items():
    print(
        f"{name:28s} "
        f"population_error={row['population_l2_error']:.6e} "
        f"density_error={row['density_frobenius_error']:.6e} "
        f"purity={row['purity']:.8f} "
        f"norm_drift={row['max_norm_drift']:.3e} "
        f"cond={row['max_condition_number']:.3e}"
    )

print("\nAcceptance:",result["acceptance"])
