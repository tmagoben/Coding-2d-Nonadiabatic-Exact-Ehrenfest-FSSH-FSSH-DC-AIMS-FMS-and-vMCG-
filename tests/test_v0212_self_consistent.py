import numpy as np

from gaussian_dynamics.synthetic_operator_provider_v21 import SyntheticLinearOperatorProviderV21,SyntheticLinearOperatorConfigV21
from gaussian_dynamics.complex_gauge_v21 import PhaseMixingGaugeV21,GaugeTransformedOperatorProviderV21,random_unitary_v21
from gaussian_dynamics.block_sparse_molecular_v21 import BlockMolecularTBFV21,BlockSparseSettingsV21
from gaussian_dynamics.block_dynamics_v21 import gauge_block_matrices_v21,gauge_mapped_coefficient_error_v21
from gaussian_dynamics.self_consistent_block_v212 import SelfConsistentBlockSettingsV212,run_self_consistent_block_dynamics_v212


def _base():
    return SyntheticLinearOperatorProviderV21(SyntheticLinearOperatorConfigV21(nstate=2,nq=2,seed=21212,mass=30.0,base_scale=.025,derivative_scale=.01))


def _basis():
    return [
        BlockMolecularTBFV21(0,np.array([-.7,.2]),np.array([.3,.05]),np.diag([1.2,1.5])),
        BlockMolecularTBFV21(1,np.array([.05,.3]),np.array([-.15,.08]),np.array([[1.7,.15],[.15,1.3]])),
        BlockMolecularTBFV21(2,np.array([.8,.15]),np.array([.2,-.04]),np.diag([1.0,1.8])),
    ]


def test_self_consistent_block_nuclear_dynamics_is_complex_gauge_equivalent():
    gauge=PhaseMixingGaugeV21(random_unitary_v21(2,21213),np.array([[.14,-.06],[-.08,.11]]),np.array([.2,-.3]))
    C0=np.array([.7+.1j,.25-.15j,.4+.2j,-.1+.3j,.3-.2j,.2+.1j])
    settings=SelfConsistentBlockSettingsV212(
        graph=BlockSparseSettingsV21(enter_score=1e-14,exit_score=1e-14,search_overlap_floor=1e-14,local_omitted_score_l2_budget=0,use_kdtree=False),
        use_dense_reference=True,corrector_iterations=3,momentum_tolerance=1e-12,
    )
    errors=[]
    for dt in (.01,.005):
        steps=int(round(.05/dt))
        base=run_self_consistent_block_dynamics_v212(_basis(),C0,_base(),dt=dt,steps=steps,settings=settings,store_every=steps)
        gp=GaugeTransformedOperatorProviderV21(_base(),gauge)
        G0,_=gauge_block_matrices_v21(gauge,_basis(),np.zeros((3,2)))
        transformed=run_self_consistent_block_dynamics_v212(_basis(),G0.conj().T@C0,gp,dt=dt,steps=steps,settings=settings,store_every=steps)
        qerr=max(np.linalg.norm(a.q-b.q) for a,b in zip(base["final_basis"],transformed["final_basis"]))
        perr=max(np.linalg.norm(a.p-b.p) for a,b in zip(base["final_basis"],transformed["final_basis"]))
        Gf,_=gauge_block_matrices_v21(gauge,base["final_basis"],np.zeros((3,2)))
        cerr=gauge_mapped_coefficient_error_v21(base["final_coefficients"],transformed["final_coefficients"],base["final_S"],Gf)
        assert qerr<1e-12
        assert perr<5e-12
        assert base["maximum_norm_drift"]<2e-10
        assert transformed["maximum_norm_drift"]<2e-10
        errors.append(cerr)
    assert errors[1]<errors[0]/3.5
    assert errors[-1]<1e-11
