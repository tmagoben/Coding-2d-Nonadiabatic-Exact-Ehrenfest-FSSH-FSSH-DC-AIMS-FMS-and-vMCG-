import numpy as np
import pytest

from gaussian_dynamics.molecular_backend import (
    MolecularGeometry,
    CartesianElectronicStructurePoint,
)
from gaussian_dynamics.pyscf_backend_v05 import PySCFSACASSCFConfig
from gaussian_dynamics.pyscf_wavefunction_overlap import CASSCFWavefunctionSnapshot
from gaussian_dynamics.pyscf_tracked_backend_v06 import (
    PySCFTrackedSACASSCFBackend,
)


def make_point(E,d01,tag):
    geom=MolecularGeometry(
        ("H","H"),
        np.array([[0.0,0.0,0.0],[0.0,0.0,1.4+0.01*tag]])
    )
    G=np.zeros((2,2,3))
    G[0,:,2]=1.0+tag
    G[1,:,2]=2.0+tag

    D=np.zeros((2,2,2,3))
    D[0,1,:,2]=d01
    D[1,0,:,2]=-d01

    point=CartesianElectronicStructurePoint(
        geometry=geom,
        energies=np.asarray(E,float),
        gradients_cart=G,
        nac_cart=D,
        masses_amu=np.array([1.0,1.0]),
        metadata={"raw_tag":tag},
    ).validate()

    roots=(
        np.array([[1.0],[0.0]]),
        np.array([[0.0],[1.0]]),
    )
    snap=CASSCFWavefunctionSnapshot(
        mol=object(),
        mo_coeff=np.eye(2),
        ci_roots=roots,
        ncore=0,
        ncas=2,
        nelecas=(1,0),
        metadata={"raw_tag":tag},
    )
    return point,snap


def test_tracked_backend_reorders_crossing_roots_and_nac():
    raw0=make_point([0.0,1.0],0.2,0)
    raw1=make_point([0.1,0.9],0.3,1)

    overlaps=[
        np.array([
            [0.05,-0.96],
            [0.94, 0.03],
        ])
    ]

    cfg=PySCFSACASSCFConfig(
        basis="fake",ncas=2,nelecas=(1,0),nstates=2
    )
    backend=PySCFTrackedSACASSCFBackend(
        cfg,
        minimum_overlap=0.5,
        minimum_score_margin=0.1,
        overlap_engine=lambda prev,curr: overlaps.pop(0),
    )

    raws=[raw0,raw1]
    backend._run_raw=lambda geom: raws.pop(0)

    dummy=raw0[0].geometry
    p0=backend.evaluate(dummy)
    p1=backend.evaluate(dummy)

    assert np.allclose(p0.energies,[0.0,1.0])
    assert np.allclose(p1.energies,[0.9,0.1])

    # Raw gradient states [0,1] are reordered to [1,0].
    assert np.allclose(p1.gradients_cart[:,0,2],[3.0,2.0])

    # raw d_10=-0.3, but state 0 receives an additional - sign:
    # (-1)*(+1)*(-0.3)=+0.3.
    assert p1.nac_cart[0,1,0,2] == pytest.approx(0.3)
    assert p1.nac_cart[1,0,0,2] == pytest.approx(-0.3)

    assert p1.metadata["permutation_tracked_to_raw"] == [1,0]
    assert p1.metadata["ambiguous"] is False


def test_tracked_backend_raises_on_ambiguous_root_identity():
    raw0=make_point([0.0,1.0],0.2,0)
    raw1=make_point([0.49,0.51],0.2,1)

    s=1/np.sqrt(2)
    overlap=np.array([[s,s],[s,-s]])

    cfg=PySCFSACASSCFConfig(
        basis="fake",ncas=2,nelecas=(1,0),nstates=2
    )
    backend=PySCFTrackedSACASSCFBackend(
        cfg,
        minimum_overlap=0.4,
        minimum_score_margin=0.01,
        ambiguity_policy="raise",
        overlap_engine=lambda prev,curr: overlap,
    )

    raws=[raw0,raw1]
    backend._run_raw=lambda geom: raws.pop(0)

    dummy=raw0[0].geometry
    backend.evaluate(dummy)

    with pytest.raises(RuntimeError,match="tracking is ambiguous"):
        backend.evaluate(dummy)


def test_reset_tracking_restarts_labels():
    raw0=make_point([0.0,1.0],0.2,0)

    cfg=PySCFSACASSCFConfig(
        basis="fake",ncas=2,nelecas=(1,0),nstates=2
    )
    backend=PySCFTrackedSACASSCFBackend(cfg)
    backend._run_raw=lambda geom: raw0

    p=backend.evaluate(raw0[0].geometry)
    assert backend.step_index==1

    backend.reset_tracking()
    assert backend.step_index==0
    assert backend.previous_snapshot is None
