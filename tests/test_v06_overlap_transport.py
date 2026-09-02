import numpy as np

from gaussian_dynamics.overlap_transport import (
    nearest_unitary,
    current_to_previous_procrustes,
    directional_nac_from_overlap,
    overlap_unitarity_defect,
)


def test_procrustes_aligns_unitary_rotated_current_basis():
    theta=0.37
    O=np.array([
        [np.cos(theta), -np.sin(theta)],
        [np.sin(theta),  np.cos(theta)],
    ],dtype=complex)

    Q,aligned,s=current_to_previous_procrustes(O)

    assert np.allclose(Q,O.conj().T,atol=1e-12)
    assert np.allclose(aligned,np.eye(2),atol=1e-12)
    assert np.allclose(s,[1.0,1.0],atol=1e-12)


def test_directional_nac_from_small_overlap_rotation():
    k=0.7
    ds=1e-5
    th=k*ds
    O=np.array([
        [np.cos(th), -np.sin(th)],
        [np.sin(th),  np.cos(th)],
    ])

    d=directional_nac_from_overlap(O,ds)
    expected=np.array([[0.0,-k],[k,0.0]])

    assert np.allclose(d,expected,atol=2e-10)
    assert overlap_unitarity_defect(O) < 1e-12


def test_nearest_unitary_removes_singular_value_distortion():
    U=np.array([[0.0,1.0],[-1.0,0.0]],dtype=complex)
    O=U @ np.diag([0.9,0.8])
    W=nearest_unitary(O)

    assert np.allclose(W.conj().T@W,np.eye(2),atol=1e-12)
