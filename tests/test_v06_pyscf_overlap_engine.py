import itertools
import numpy as np

import gaussian_dynamics.pyscf_wavefunction_overlap as pwo
from gaussian_dynamics.pyscf_wavefunction_overlap import (
    CASSCFWavefunctionSnapshot,
    embed_active_ci_with_doubly_occupied_core,
    casscf_state_overlap_matrix,
)


class FakeCistring:
    @staticmethod
    def make_strings(orb_list, nelec):
        orb_list=list(orb_list)
        out=[]
        for occ in itertools.combinations(orb_list,nelec):
            mask=0
            for i in occ:
                mask |= (1 << i)
            out.append(mask)
        return np.asarray(out,dtype=np.int64)

    @classmethod
    def strs2addr(cls,norb,nelec,strings):
        base=cls.make_strings(range(norb),nelec)
        lookup={int(s):i for i,s in enumerate(base)}
        return np.asarray([lookup[int(s)] for s in strings],dtype=int)


def occs(mask,norb):
    return [i for i in range(norb) if int(mask) & (1 << i)]


class FakeAddons:
    @staticmethod
    def overlap(bra,ket,norb,nelec,s=None):
        na,nb=nelec
        strsa=FakeCistring.make_strings(range(norb),na)
        strsb=FakeCistring.make_strings(range(norb),nb)

        total=0.0+0.0j
        for ia,sa in enumerate(strsa):
            oa=occs(sa,norb)
            for ib,sb in enumerate(strsb):
                ob=occs(sb,norb)
                cb=np.conj(bra[ia,ib])
                if abs(cb)==0:
                    continue
                for ja,ta in enumerate(strsa):
                    pa=occs(ta,norb)
                    det_a=np.linalg.det(s[np.ix_(oa,pa)]) if na else 1.0
                    for jb,tb in enumerate(strsb):
                        ck=ket[ja,jb]
                        if abs(ck)==0:
                            continue
                        pb=occs(tb,norb)
                        det_b=np.linalg.det(s[np.ix_(ob,pb)]) if nb else 1.0
                        total += cb*ck*det_a*det_b
        return total


class FakeFCI:
    addons=FakeAddons()


class FakeGTO:
    @staticmethod
    def intor_cross(name,mol1,mol2):
        assert name=="int1e_ovlp_sph"
        return np.asarray(mol1.cross_overlap)


class FakeMol:
    def __init__(self,cross_overlap):
        self.cross_overlap=np.asarray(cross_overlap)


def fake_imports():
    return FakeGTO,FakeFCI,FakeCistring


def test_active_ci_embedding_adds_doubly_occupied_core():
    ci=np.array([[1.0]])
    full=embed_active_ci_with_doubly_occupied_core(
        ci,
        ncore=1,
        ncas=1,
        nelecas=(1,0),
        cistring=FakeCistring,
    )

    # alpha has core + active occupied; beta has only core occupied.
    # There is only one allowed determinant for this tiny space.
    assert full.shape==(1,2)
    assert np.count_nonzero(full)==1
    assert np.max(np.abs(full))==1.0


def test_full_core_active_overlap_includes_cross_subspace_mixing(monkeypatch):
    monkeypatch.setattr(pwo,"_imports",fake_imports)

    theta=0.31
    c=np.cos(theta)
    s=np.sin(theta)

    Cprev=np.eye(2)
    Ccurr=np.array([[c,-s],[s,c]])

    # AO cross overlap is identity, so MO cross overlap is Ccurr.
    mol_prev=FakeMol(np.eye(2))
    mol_curr=FakeMol(np.eye(2))

    prev=CASSCFWavefunctionSnapshot(
        mol=mol_prev,
        mo_coeff=Cprev,
        ci_roots=(np.array([[1.0]]),),
        ncore=1,
        ncas=1,
        nelecas=(1,0),
        metadata={},
    )
    curr=CASSCFWavefunctionSnapshot(
        mol=mol_curr,
        mo_coeff=Ccurr,
        ci_roots=(np.array([[1.0]]),),
        ncore=1,
        ncas=1,
        nelecas=(1,0),
        metadata={},
    )

    O=casscf_state_overlap_matrix(prev,curr)

    # Alpha occupies both correlated orbitals -> det(rotation)=1.
    # Beta occupies only the core orbital -> overlap cos(theta).
    # The exact core+active many-electron overlap is therefore cos(theta).
    np.testing.assert_allclose(O[0,0],c,atol=1e-12)


def test_many_electron_overlap_detects_root_swap_and_sign(monkeypatch):
    monkeypatch.setattr(pwo,"_imports",fake_imports)

    mol_prev=FakeMol(np.eye(2))
    mol_curr=FakeMol(np.eye(2))

    e0=np.array([[1.0],[0.0]])
    e1=np.array([[0.0],[1.0]])

    prev=CASSCFWavefunctionSnapshot(
        mol=mol_prev,
        mo_coeff=np.eye(2),
        ci_roots=(e0,e1),
        ncore=0,
        ncas=2,
        nelecas=(1,0),
        metadata={},
    )
    curr=CASSCFWavefunctionSnapshot(
        mol=mol_curr,
        mo_coeff=np.eye(2),
        ci_roots=(e1,-e0),
        ncore=0,
        ncas=2,
        nelecas=(1,0),
        metadata={},
    )

    O=casscf_state_overlap_matrix(prev,curr)
    expected=np.array([[0.0,-1.0],[1.0,0.0]])

    assert np.allclose(O,expected,atol=1e-12)
