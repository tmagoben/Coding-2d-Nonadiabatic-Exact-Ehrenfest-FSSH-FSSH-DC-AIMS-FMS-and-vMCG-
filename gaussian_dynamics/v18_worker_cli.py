import argparse
from pathlib import Path

from .convergence_worker_v18 import (
    run_coordinate_worker_v18,
)
from .campaign_io import save_campaign_json


def main():
    parser=argparse.ArgumentParser(
        description="Run one isolated v0.18 convergence coordinate."
    )
    parser.add_argument("--output",required=True)
    parser.add_argument("--dt",type=float,required=True)
    parser.add_argument("--max-basis",type=int,required=True)
    parser.add_argument("--edge-budget",type=float,required=True)
    parser.add_argument("--enrich-threshold",type=float,required=True)
    parser.add_argument("--final-time",type=float,default=0.60)
    parser.add_argument("--trajectory",action="store_true")
    parser.add_argument(
        "--trajectory-store-interval",
        type=float,default=0.10,
    )
    args=parser.parse_args()

    result=run_coordinate_worker_v18(
        {
            "dt":args.dt,
            "max_basis":args.max_basis,
            "local_score_budget":args.edge_budget,
            "enrich_threshold":args.enrich_threshold,
        },
        final_time=args.final_time,
        trajectory=args.trajectory,
        trajectory_store_interval=
            args.trajectory_store_interval,
    )
    path=save_campaign_json(
        Path(args.output),result
    )

    row=result["result"]
    print("saved:",path)
    print(
        "projected fidelity:",
        row["projected_wavefunction_fidelity"],
    )
    print("basis:",row["basis_size"])


if __name__=="__main__":
    main()
