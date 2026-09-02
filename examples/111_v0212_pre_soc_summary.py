from pathlib import Path
import json
root=Path(__file__).resolve().parents[1]
data=json.loads((root/'results/v0212_pre_soc_hardening_campaign.json').read_text())
print('v0.21.2 pre-SOC integration hardening')
print('acceptance:',data['acceptance']['passed'])
print('checks:',data['acceptance']['checks'])
print('PySCF runtime validated:',data['pyscf']['runtime_validated'])
