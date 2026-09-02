import numpy as np

from gaussian_dynamics.block_sparse_molecular_v21 import BlockMolecularTBFV21
from gaussian_dynamics.block_basis_lifecycle_v212 import insert_zero_block_v212,prune_block_projected_v212


def test_zero_block_birth_preserves_state_and_projection_prune_matches_schur_loss():
    s=2
    A=np.eye(2)
    basis=[
        BlockMolecularTBFV21(0,np.array([-0.5,0.0]),np.zeros(2),A),
        BlockMolecularTBFV21(1,np.array([0.5,0.0]),np.zeros(2),A),
    ]
    C=np.array([0.8+0.1j,0.2-0.3j,-0.1+0.2j,0.4+0.05j])
    new=BlockMolecularTBFV21(2,np.array([1.2,0.0]),np.zeros(2),1.4*A)
    born,Cborn=insert_zero_block_v212(basis,C,new,s)
    assert len(born)==3
    assert np.allclose(Cborn[:4],C)
    assert np.allclose(Cborn[4:],0.0)

    # Positive-definite block metric with nonzero retained/deleted coupling.
    X=np.array([
        [1.0,0.0,0.10+0.03j,0.0,0.04,0.0],
        [0.0,1.0,0.0,0.08-0.02j,0.0,0.03],
        [0.10-0.03j,0.0,1.0,0.0,0.12,0.0],
        [0.0,0.08+0.02j,0.0,1.0,0.0,0.10],
        [0.04,0.0,0.12,0.0,1.0,0.0],
        [0.0,0.03,0.0,0.10,0.0,1.0],
    ],dtype=complex)
    S=X.conj().T@X
    result=prune_block_projected_v212(born,Cborn,S,s,2)
    assert result.removed_uid==2
    assert result.projection_loss<1e-20  # deleted block coefficient was exactly zero
    assert np.allclose(result.coefficients,C,atol=1e-12)


def test_self_consistent_runner_accepts_zero_block_birth_at_step_boundary():
    from gaussian_dynamics.synthetic_operator_provider_v21 import SyntheticLinearOperatorConfigV21,SyntheticLinearOperatorProviderV21
    from gaussian_dynamics.block_sparse_molecular_v21 import BlockSparseSettingsV21
    from gaussian_dynamics.self_consistent_block_v212 import SelfConsistentBlockSettingsV212,run_self_consistent_block_dynamics_v212

    provider=SyntheticLinearOperatorProviderV21(SyntheticLinearOperatorConfigV21(nstate=2,nq=2,seed=21250,mass=20.0))
    A=np.eye(2)
    basis=[
        BlockMolecularTBFV21(0,np.array([-.4,.1]),np.array([.1,0]),A),
        BlockMolecularTBFV21(1,np.array([.4,.1]),np.array([-.1,0]),1.2*A),
    ]
    C0=np.array([.8+.1j,.2-.1j,.3+.1j,-.1+.2j])
    newborn=BlockMolecularTBFV21(9,np.array([1.8,.1]),np.array([0.,0.]),1.5*A)
    def policy(step,basis,C,S):
        return {"insert":newborn} if step==1 else None
    settings=SelfConsistentBlockSettingsV212(
        graph=BlockSparseSettingsV21(enter_score=1e-14,exit_score=1e-14,search_overlap_floor=1e-14,local_omitted_score_l2_budget=0,use_kdtree=False),
        use_dense_reference=True,corrector_iterations=2,
    )
    out=run_self_consistent_block_dynamics_v212(basis,C0,provider,dt=.001,steps=1,settings=settings,store_every=1,adaptation_policy=policy)
    assert out["adaptation_events"][0]["kind"]=="insert"
    assert len(out["final_basis"])==3
    assert out["final_coefficients"].shape==(6,)
    assert np.allclose(out["final_coefficients"][-2:],0.0)
    assert out["maximum_norm_drift"]<1e-8
