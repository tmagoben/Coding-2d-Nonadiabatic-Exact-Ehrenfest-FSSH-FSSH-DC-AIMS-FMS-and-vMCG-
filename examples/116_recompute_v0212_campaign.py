from pathlib import Path
from gaussian_dynamics import run_v0212_release_benchmark
from gaussian_dynamics.campaign_io import save_campaign_json
root=Path(__file__).resolve().parents[1]
out=run_v0212_release_benchmark()
path=save_campaign_json(root/'results/v0212_pre_soc_hardening_campaign_recomputed.json',out)
print('Saved:',path)
print('Acceptance:',out['acceptance'])
