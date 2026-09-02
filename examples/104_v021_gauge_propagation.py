from pathlib import Path
import json
root=Path(__file__).resolve().parents[1]; d=json.loads((root/'results'/'v021_complex_block_framework_campaign.json').read_text())
for r in d['gauge_propagation']['rows']: print(r)
print('orders',d['gauge_propagation']['observed_orders'])
