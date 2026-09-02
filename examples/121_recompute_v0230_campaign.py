from pathlib import Path

from gaussian_dynamics import (
    build_v230_reference_replay,
    build_v230_doublet_reference_replay,
    run_v0230_release_benchmark,
)
from gaussian_dynamics.campaign_io import save_campaign_json


root = Path(__file__).resolve().parents[1]
replay = build_v230_reference_replay(
    root / "results/v0230_reference_replay/even_singlet_triplet",
    overwrite=True,
)
doublet_replay = build_v230_doublet_reference_replay(
    root / "results/v0230_reference_replay/odd_doublet",
    overwrite=True,
)
output = run_v0230_release_benchmark()
path = save_campaign_json(
    root / "results/v0230_molecular_soc_admission_campaign.json", output
)
print("Replay fingerprint:", replay.dataset_fingerprint)
print("Doublet replay fingerprint:", doublet_replay.dataset_fingerprint)
print("Saved:", path)
print("Acceptance:", output["acceptance"])
print("Claims:", output["claims"])
