from pathlib import Path
import json
root=Path(__file__).resolve().parents[1]
data=json.loads((root/'results/v0212_pre_soc_hardening_campaign.json').read_text())
for row in data['self_consistent_block_dynamics']['rows']:
    print(row)
print('orders:',data['self_consistent_block_dynamics']['observed_orders'])
