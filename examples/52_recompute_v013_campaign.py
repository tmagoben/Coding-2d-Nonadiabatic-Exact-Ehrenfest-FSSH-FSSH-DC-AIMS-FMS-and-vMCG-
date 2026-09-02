from pathlib import Path

from gaussian_dynamics import (
    run_v013_release_benchmark,
    save_campaign_json,
)

root=Path(__file__).resolve().parents[1]

result=run_v013_release_benchmark(
    repository_root=root,
)

path=save_campaign_json(
    root/"results"/"v013_residual_driven_campaign_recomputed.json",
    result,
)

print("Recomputed v0.13 residual-driven campaign.")
print("Saved:",path)
print("Acceptance:",result["acceptance"])
print(
    "Reference target-density error:",
    result["reference"]["target_density_error"],
)
print(
    "Projected-state dynamics error:",
    result["reference"]["projected_dynamics_density_error"],
)
