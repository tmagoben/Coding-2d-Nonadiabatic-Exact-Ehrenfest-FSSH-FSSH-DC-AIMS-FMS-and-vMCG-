from pathlib import Path
import json
root=Path(__file__).resolve().parents[1]; d=json.loads((root/'results'/'v021_complex_block_framework_campaign.json').read_text())
print('threshold'); [print(r) for r in d['block_sparse_convergence']['threshold_rows']]
print('budget'); [print(r) for r in d['block_sparse_convergence']['budget_rows']]
