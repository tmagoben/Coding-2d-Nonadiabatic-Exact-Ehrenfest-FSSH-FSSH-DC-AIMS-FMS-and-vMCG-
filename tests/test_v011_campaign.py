import json
import numpy as np

from gaussian_dynamics.benchmark_campaign import CIPassageConfig
from gaussian_dynamics.v11_benchmark import run_v011_release_benchmark
from gaussian_dynamics.campaign_io import save_campaign_json


def test_tiny_v011_release_campaign_runs_without_ablations(tmp_path):
    config=CIPassageConfig(
        q0=(0.55,0.45),
        p0=(0.6,0.8),
        A_diag=(1.2,1.2),
        mass=20.0,
        final_time=0.01,
        half_width=3.0,
    )

    result=run_v011_release_benchmark(
        config,
        include_ablations=False,
    )

    assert len(result["basis_ladder"])==5
    assert "v11_reference" in result
    assert "acceptance" in result


def test_campaign_json_serializes_complex_density(tmp_path):
    payload={
        "rho":np.array([
            [0.7+0j,0.1+0.2j],
            [0.1-0.2j,0.3+0j],
        ]),
    }
    path=save_campaign_json(tmp_path/"campaign.json",payload)
    loaded=json.loads(path.read_text())

    assert loaded["rho"][0][1]==[0.1,0.2]
