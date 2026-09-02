from pathlib import Path

from gaussian_dynamics import (
    build_v231_admission_bundle,
    run_v0231_release_benchmark,
)
from gaussian_dynamics.campaign_io import save_campaign_json


root = Path(__file__).resolve().parents[1]
bundle_root = root / "results/v0231_admission_bundles"
even = build_v231_admission_bundle(
    bundle_root / "even_singlet_triplet", overwrite=True
)
odd = build_v231_admission_bundle(
    bundle_root / "odd_doublet", odd=True, overwrite=True
)
output = run_v0231_release_benchmark()
path = save_campaign_json(
    root / "results/v0231_raw_evidence_admission_campaign.json", output
)
print("Even replay fingerprint:", even["dataset"].dataset_fingerprint)
print("Even dossier fingerprint:", even["dossier"].fingerprint())
print("Odd replay fingerprint:", odd["dataset"].dataset_fingerprint)
print("Odd dossier fingerprint:", odd["dossier"].fingerprint())
print("Saved:", path)
print("Acceptance:", output["acceptance"])
print("Claims:", output["claims"])
