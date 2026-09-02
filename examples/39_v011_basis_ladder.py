from pathlib import Path
import json

root=Path(__file__).resolve().parents[1]
result=json.loads(
    (root/"results"/"v011_basis_completeness_campaign.json").read_text()
)

print("v0.11 saved basis ladder")
print("------------------------")
print("exact populations:",result["exact"]["populations"])

for row in result["basis_ladder"]:
    print(
        f"Nmax={row['max_basis']:2d} "
        f"Nfinal={row['basis_size']:2d} "
        f"population_error={row['population_l2_error']:.6e} "
        f"density_error={row['density_frobenius_error']:.6e} "
        f"purity={row['purity']:.8f} "
        f"norm_drift={row['max_norm_drift']:.3e} "
        f"cond={row['max_condition_number']:.3e}"
    )

print(
    "\nThis example reads the release campaign instead of recomputing several "
    "large dynamic gauge graphs. Use example 41 to regenerate the base campaign."
)
