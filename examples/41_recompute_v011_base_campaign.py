from pathlib import Path

from gaussian_dynamics import (
    run_v011_release_benchmark,
    save_campaign_json,
)

root=Path(__file__).resolve().parents[1]

result=run_v011_release_benchmark(
    include_ablations=False,
)

path=save_campaign_json(
    root/"results"/"v011_basis_completeness_campaign_recomputed_base.json",
    result,
)

print("Recomputed the exact reference, v0.10 baseline, and v0.11 basis ladder.")
print("Saved:",path)
print("Reference population error:",
      result["v11_reference"]["population_l2_error"])
print("Reference full-density error:",
      result["v11_reference"]["density_frobenius_error"])
print("Acceptance:",result["acceptance"])
print(
    "\nThe release ablations are intentionally run in separate processes to avoid "
    "retaining several large dynamic gauge graphs simultaneously."
)
