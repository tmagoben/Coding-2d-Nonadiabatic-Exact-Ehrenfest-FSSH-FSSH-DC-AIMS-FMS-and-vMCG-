import numpy as np

from gaussian_dynamics.dynamic_graph_aims import DynamicGraphTBF
from gaussian_dynamics.managed_graph_aims import run_managed_graph_aims


def test_managed_graph_aims_spawns_and_preserves_norm():
    b=DynamicGraphTBF(
        0,1,np.array([0.55,0.45]),np.array([0.6,0.8]),1.2*np.eye(2),('seed',0)
    )
    out=run_managed_graph_aims(
        [b],[1+0j],dt=2e-4,steps=30,spa_order=1,
        spawn_action_threshold=1e-6,max_basis=2,store_every=5
    )
    spawn=[e for e in out['events'] if e['kind']=='spawn']
    assert len(spawn)==1
    norms=np.array([r['norm'] for r in out['records']])
    assert np.max(np.abs(norms-1.0)) < 5e-5
    assert abs(out['final_coefficients'][1]) > 1e-10
    assert all(np.isfinite(r['spa1_relative_correction']) for r in out['records'])

def test_managed_graph_aims_prunes_initial_redundancy_before_propagation():
    A=1.2*np.eye(2)
    b1=DynamicGraphTBF(0,1,np.array([0.55,0.45]),np.array([0.6,0.8]),A,('seed',0))
    b2=DynamicGraphTBF(1,1,np.array([0.550000001,0.45]),np.array([0.6,0.8]),A,('seed',1))

    out=run_managed_graph_aims(
        [b1,b2],[0.5+0j,0.5+0j],
        dt=2e-4,steps=2,spa_order=0,
        spawn_action_threshold=1.0,max_basis=2,
        condition_limit=1e6,eigenvalue_floor=1e-8,
        max_pruning_loss=1e-5,store_every=1,
    )

    prunes=[e for e in out['events'] if e['kind']=='prune']
    assert len(prunes)==1
    assert prunes[0]['step']==0
    assert len(out['final_basis'])==1
    assert max(r['condition_number'] for r in out['records']) < 2.0
