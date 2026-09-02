from pathlib import Path
from gaussian_dynamics import run_v021_release_benchmark
from gaussian_dynamics.campaign_io import save_campaign_json
root=Path(__file__).resolve().parents[1]; out=run_v021_release_benchmark(); path=save_campaign_json(root/'results'/'v021_complex_block_framework_campaign_recomputed.json',out); print(path); print(out['acceptance'])
