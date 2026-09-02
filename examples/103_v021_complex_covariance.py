from pathlib import Path
import json
root=Path(__file__).resolve().parents[1]
data=json.loads((root/'results'/'v021_complex_block_framework_campaign.json').read_text())
print(data['point_covariance']); print(data['block_covariance'])
