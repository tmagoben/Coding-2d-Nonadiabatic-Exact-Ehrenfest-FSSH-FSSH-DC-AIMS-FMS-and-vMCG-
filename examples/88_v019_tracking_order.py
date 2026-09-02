from pathlib import Path
import json

root=Path(__file__).resolve().parents[1]
data=json.loads(
    (root/"results"/"v019_molecular_direct_dynamics_campaign.json").read_text()
)
p=data["provider_scan"]

print("v0.19 nearest-anchor state tracking")
print("-----------------------------------")
print("after reference seed, shuffled query order:")
print("max energy error:",
      p["shuffled_after_seed_max_energy_error"])
print("max NAC error:",
      p["shuffled_after_seed_max_nac_error"])
print("ambiguities:",
      p["shuffled_after_seed_diagnostics"]["tracking_ambiguities"])
print("cache size:",
      p["scrambled_provider_diagnostics"]["cache_size"])
print("cache hits:",
      p["scrambled_provider_diagnostics"]["cache_hits"])
