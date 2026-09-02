from pathlib import Path
import json

root=Path(__file__).resolve().parents[1]
data=json.loads(
    (root/"results"/"v013_residual_driven_campaign.json").read_text()
)

pure=data["pure_residual_greedy"]
screened=data["density_screened_residual"]

print("v0.13 residual-driven initial basis")
print("-----------------------------------")
print("pure residual greedy:")
print("  basis size:",pure["basis_size"])
print("  projection fidelity:",pure["projection_fidelity"])
print("  relative residual:",pure["relative_residual"])
print("  initial density error:",pure["initial_density_error"])
print("  condition number:",pure["condition_number"])

print("\ndensity-screened residual reference:")
print("  basis size:",screened["basis_size"])
print("  projection fidelity:",screened["projection_fidelity"])
print("  relative residual:",screened["relative_residual"])
print("  initial density error:",screened["initial_density_error"])
print("  condition number:",screened["condition_number"])

print("\nSelected Gaussian sequence:")
for i,label in enumerate(screened["selected_labels"],start=2):
    print(f"  N={i:2d}: {label}")
