from pathlib import Path

from gaussian_dynamics import run_v0214_release_benchmark
from gaussian_dynamics.campaign_io import save_campaign_json


root = Path(__file__).resolve().parents[1]
output = run_v0214_release_benchmark()
path = save_campaign_json(
    root / "results/v0214_differential_restart_campaign.json", output
)
print("Saved:", path)
print("Acceptance:", output["acceptance"])
