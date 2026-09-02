from pathlib import Path
import json
import subprocess
import sys
import tempfile

from gaussian_dynamics.v18_benchmark import (
    sampled_audit_scaling_v18,
    assemble_v018_campaign_from_partials,
)
from gaussian_dynamics.campaign_io import (
    save_campaign_json,
)


root=Path(__file__).resolve().parents[1]

jobs={
    "canonical":(0.005,12,0.02,0.020,0.60,True,0.10),
    "dtc":      (0.010,12,0.02,0.020,0.60,False,0.10),
    "dtf":      (0.0025,12,0.02,0.020,0.60,False,0.10),
    "b10":      (0.005,10,0.02,0.020,0.60,False,0.10),
    "b11":      (0.005,11,0.02,0.020,0.60,False,0.10),
    "e08":      (0.005,12,0.08,0.020,0.60,False,0.10),
    "e0":       (0.005,12,0.00,0.020,0.60,False,0.10),
    "th32":     (0.005,12,0.02,0.032,0.60,False,0.10),
    "th28":     (0.005,12,0.02,0.028,0.60,False,0.10),
    "coarse":   (0.010,10,0.08,0.040,0.60,False,0.10),
    "fine":     (0.0025,12,0.00,0.020,0.60,False,0.10),
    "long":     (0.005,12,0.02,0.020,1.20,True,0.20),
}


def run_job(name,spec,outfile):
    dt,max_basis,budget,threshold,final_time,trajectory,interval=spec
    cmd=[
        sys.executable,
        "-m","gaussian_dynamics.v18_worker_cli",
        "--output",str(outfile),
        "--dt",str(dt),
        "--max-basis",str(max_basis),
        "--edge-budget",str(budget),
        "--enrich-threshold",str(threshold),
        "--final-time",str(final_time),
        "--trajectory-store-interval",str(interval),
    ]
    if trajectory:
        cmd.append("--trajectory")

    # A fresh interpreter is deliberate. Retry once if an external BLAS/allocator
    # stall makes a coordinate exceed the generous timeout.
    last=None
    for attempt in range(2):
        try:
            subprocess.run(
                cmd,
                cwd=root,
                check=True,
                timeout=120,
            )
            return
        except (subprocess.TimeoutExpired,subprocess.CalledProcessError) as exc:
            last=exc
            print(f"{name}: attempt {attempt+1} failed: {exc}")
    raise RuntimeError(f"{name} failed after retry") from last


with tempfile.TemporaryDirectory(
    prefix="v018_campaign_"
) as tmp:
    tmp=Path(tmp)
    parts={}

    for name,spec in jobs.items():
        print("running:",name)
        path=tmp/f"{name}.json"
        run_job(name,spec,path)
        parts[name]=json.loads(
            path.read_text(encoding="utf-8")
        )

    sampled=sampled_audit_scaling_v18()

    campaign=assemble_v018_campaign_from_partials(
        parts,
        sampled,
        repository_root=root,
    )

    path=save_campaign_json(
        root/"results"/"v018_convergence_complete_campaign_recomputed.json",
        campaign,
    )

print("saved:",path)
print("acceptance:",campaign["acceptance"])
