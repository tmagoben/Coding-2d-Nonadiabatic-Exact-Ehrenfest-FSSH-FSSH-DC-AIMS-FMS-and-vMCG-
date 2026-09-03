"""Recompute the v0.28.0 moving-frame development evidence."""
import json
from pathlib import Path
from gaussian_dynamics.moving_frame_evidence_v280 import run_moving_frame_evidence_v280

def main():
    evidence=run_moving_frame_evidence_v280(); payload=evidence.as_dict(); output=Path('results/v0280_moving_frame_evidence.json'); output.parent.mkdir(parents=True,exist_ok=True); output.write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'passed':evidence.passed,'check_count':evidence.check_count,'fingerprint':evidence.fingerprint(),'output':str(output)},indent=2))
if __name__=='__main__': main()
