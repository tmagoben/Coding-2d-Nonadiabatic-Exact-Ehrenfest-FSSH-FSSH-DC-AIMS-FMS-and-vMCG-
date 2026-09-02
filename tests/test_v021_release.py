from gaussian_dynamics import run_v021_release_benchmark

def test_v021_release_benchmark_passes():
    out=run_v021_release_benchmark(); assert out['acceptance']['passed']; assert out['nstate_scaling']['rows'][-1]['nstate']==8; assert out['dynamic_topology']['total_entered_edges']>=10; assert out['dynamic_topology']['total_exited_edges']>=5
