import numpy as np
import pytest

from gaussian_dynamics.molecular_backend import MolecularGeometry
from gaussian_dynamics.pyscf_backend_v05 import (
    PySCFSACASSCFConfig,
    PySCFSACASSCFBackend,
)
from gaussian_dynamics.pyscf_provider import PySCFStateAveragedCASSCFProvider
from gaussian_dynamics.pyscf_tracked_backend_v06 import (
    PySCFTrackedSACASSCFBackend,
)


class FakePyscf:
    __version__="2.13.1"


class FakeMol:
    def __init__(self, kwargs, log):
        self.kwargs=kwargs
        self.log=log
        self.natm=len(kwargs["atom"])

    def atom_mass_list(self,isotope_avg=False):
        self.log.append(("atom_mass_list", isotope_avg))
        return np.array([7.0,1.0])


class FakeGto:
    def __init__(self,log):
        self.log=log

    def M(self,**kwargs):
        self.log.append(("gto.M", kwargs))
        return FakeMol(kwargs,self.log)


class FakeMF:
    def __init__(self,mol,log,kind):
        self.mol=mol
        self.log=log
        self.kind=kind
        self.conv_tol=None
        self.max_cycle=None
        self.converged=True

    def kernel(self):
        self.log.append(("scf.kernel",self.kind,self.conv_tol,self.max_cycle))


class FakeScf:
    def __init__(self,log):
        self.log=log

    def RHF(self,mol):
        self.log.append(("RHF",))
        return FakeMF(mol,self.log,"RHF")

    def ROHF(self,mol):
        self.log.append(("ROHF",))
        return FakeMF(mol,self.log,"ROHF")


class FakeGrad:
    def __init__(self,state,log,natom):
        self.state=state
        self.log=log
        self.natom=natom

    def kernel(self):
        self.log.append(("grad.kernel",self.state))
        return np.full((self.natom,3),0.01*(self.state+1))


class FakeNAC:
    def __init__(self,log,natom):
        self.log=log
        self.natom=natom

    def kernel(self,state,use_etfs=False,mult_ediff=False):
        self.log.append(("nac.kernel",state,use_etfs,mult_ediff))
        ket,bra=state
        base=(10*bra+ket+1)*0.001
        if mult_ediff:
            base*=5.0
        return np.full((self.natom,3),base)


class FakeMC:
    def __init__(self,mf,ncas,nelecas,log):
        self.mf=mf
        self.ncas=ncas
        self.nelecas=nelecas
        self.log=log
        self.conv_tol=None
        self.conv_tol_grad=None
        self.max_cycle_macro=None
        self.converged=True
        self.e_states=np.array([-7.8,-7.7])
        self.mo_coeff=np.eye(4)
        self.ncore=0
        if isinstance(nelecas, tuple):
            self.nelecas=nelecas
        else:
            self.nelecas=((int(nelecas)+1)//2, int(nelecas)//2)
        self.ci=(np.eye(2),np.fliplr(np.eye(2)))

    def state_average_(self,weights):
        self.log.append(("state_average_",tuple(np.asarray(weights))))
        return self

    def state_average(self,weights):
        self.log.append(("state_average",tuple(np.asarray(weights))))
        return self

    def kernel(self,*args):
        self.log.append(("mc.kernel",len(args)))

    def nuc_grad_method(self,state=None):
        self.log.append(("nuc_grad_method",state))
        return FakeGrad(state,self.log,self.mf.mol.natm)

    def nac_method(self):
        self.log.append(("nac_method",))
        return FakeNAC(self.log,self.mf.mol.natm)


class FakeMcscf:
    def __init__(self,log):
        self.log=log

    def CASSCF(self,mf,ncas,nelecas):
        self.log.append(("CASSCF",ncas,nelecas))
        return FakeMC(mf,ncas,nelecas,self.log)


def make_backend(log,**kwargs):
    cfg=PySCFSACASSCFConfig(
        basis="sto-3g",
        ncas=2,
        nelecas=2,
        nstates=2,
        use_etfs=True,
        compute_scaled_nac=True,
        **kwargs,
    )
    backend=PySCFSACASSCFBackend(cfg)
    backend._imports=lambda: (
        FakePyscf,
        FakeGto(log),
        FakeScf(log),
        FakeMcscf(log),
    )
    return backend


def test_explicit_pyscf_backend_call_contract():
    log=[]
    backend=make_backend(log)

    geom=MolecularGeometry(
        ("Li","H"),
        np.array([[0.0,0.0,0.0],[0.0,0.0,3.0]])
    )

    point=backend.evaluate(geom)

    gto_call=[x for x in log if x[0]=="gto.M"][0]
    kwargs=gto_call[1]

    assert kwargs["unit"]=="Bohr"
    assert kwargs["basis"]=="sto-3g"
    assert ("RHF",) in log
    assert any(x[0]=="state_average_" and np.allclose(x[1],[0.5,0.5]) for x in log)
    assert ("nuc_grad_method",0) in log
    assert ("nuc_grad_method",1) in log

    # v0.23.2 freezes the empirically certified PySCF 2.13.1 mapping:
    # internal d[0,1]=<0|d 1> is requested with state=(0,1).
    assert ("nac.kernel",(0,1),True,False) in log
    assert ("nac.kernel",(0,1),True,True) in log

    assert np.allclose(point.nac_cart[1,0],-point.nac_cart[0,1])
    assert np.allclose(point.scaled_nac_cart[1,0],point.scaled_nac_cart[0,1])
    assert point.metadata["pyscf_nac_state_tuple"]=="(ket,bra)"
    assert point.metadata["pyscf_request_for_internal_dij"]=="state=(i,j)"
    assert "central differences" in point.metadata["pyscf_nac_empirical_mapping"]
    assert point.metadata["dynamics_mult_ediff"] is False


def test_rohf_reference_is_explicitly_selectable():
    log=[]
    backend=make_backend(log,scf_reference="ROHF")

    geom=MolecularGeometry(
        ("Li","H"),
        np.array([[0.0,0.0,0.0],[0.0,0.0,3.0]])
    )
    backend.evaluate(geom)

    assert ("ROHF",) in log
    assert ("RHF",) not in log


def test_tracked_backend_uses_empirical_v232_nac_mapping():
    log=[]
    cfg=PySCFSACASSCFConfig(
        basis="sto-3g",
        ncas=2,
        nelecas=2,
        nstates=2,
        use_etfs=False,
        compute_scaled_nac=True,
    )
    backend=PySCFTrackedSACASSCFBackend(cfg)
    backend._imports=lambda: (
        FakePyscf,
        FakeGto(log),
        FakeScf(log),
        FakeMcscf(log),
    )
    geom=MolecularGeometry(
        ("Li","H"),
        np.array([[0.0,0.0,0.0],[0.0,0.0,3.0]])
    )

    point,_=backend.evaluate_raw_with_snapshot(geom)

    assert ("nac.kernel",(0,1),False,False) in log
    assert ("nac.kernel",(0,1),False,True) in log
    assert point.metadata["pyscf_request_for_internal_dij"]=="state=(i,j)"


def test_one_dimensional_provider_uses_empirical_v232_nac_mapping():
    log=[]

    def builder(q):
        atoms=[("Li",(0.0,0.0,0.0)),("H",(0.0,0.0,q))]
        tangent=np.array([[0.0,0.0,0.0],[0.0,0.0,1.0]])
        return atoms,tangent

    provider=PySCFStateAveragedCASSCFProvider(
        builder,
        basis="sto-3g",
        ncas=2,
        nelecas=2,
        nstates=2,
        use_etfs=False,
    )
    provider._imports=lambda: (FakeGto(log),FakeScf(log),FakeMcscf(log))

    point=provider.evaluate(3.0)

    assert ("nac.kernel",(0,1),False,False) in log
    assert point.metadata["pyscf_request_for_internal_dij"]=="state=(i,j)"


def test_missing_real_pyscf_has_clear_error_when_not_installed():
    cfg=PySCFSACASSCFConfig(
        basis="sto-3g",ncas=2,nelecas=2,nstates=2
    )
    backend=PySCFSACASSCFBackend(cfg)

    try:
        import pyscf  # noqa
    except ImportError:
        geom=MolecularGeometry(
            ("Li","H"),
            np.array([[0.0,0.0,0.0],[0.0,0.0,3.0]])
        )
        with pytest.raises(ImportError,match="PySCF is required"):
            backend.evaluate(geom)
