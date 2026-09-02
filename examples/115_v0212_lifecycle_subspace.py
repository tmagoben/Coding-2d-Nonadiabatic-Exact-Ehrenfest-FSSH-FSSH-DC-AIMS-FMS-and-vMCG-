from pathlib import Path
import json
root=Path(__file__).resolve().parents[1]
data=json.loads((root/'results/v0212_pre_soc_hardening_campaign.json').read_text())
print('lifecycle:',data['adaptive_block_lifecycle'])
print('subspace checks:',data['subspace_provider']['subspace_checks'])
print('minimum singular value:',data['subspace_provider']['minimum_seen_singular_value'])
