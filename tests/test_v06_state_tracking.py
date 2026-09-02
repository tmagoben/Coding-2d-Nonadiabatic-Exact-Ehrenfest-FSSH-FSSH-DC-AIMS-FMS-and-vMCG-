import numpy as np
import pytest

from gaussian_dynamics.state_tracking import (
    maximum_overlap_assignment,
    transform_state_properties,
    energy_degeneracy_clusters,
    subspace_overlap_singular_values,
)


def test_maximum_overlap_assignment_recovers_swap_and_sign():
    O=np.array([
        [0.05, -0.96],
        [0.94,  0.03],
    ])

    result=maximum_overlap_assignment(
        O,
        minimum_overlap=0.5,
        minimum_score_margin=0.1,
        real_gauge=True,
    )

    assert np.array_equal(result.permutation,[1,0])
    assert np.array_equal(np.real(result.phase_factors),[-1.0,1.0])
    assert not result.ambiguous
    assert np.all(np.real(result.assigned_overlaps) > 0.0)


def test_state_property_transform_reorders_and_gauge_transforms_nac():
    O=np.array([
        [0.05, -0.96],
        [0.94,  0.03],
    ])
    result=maximum_overlap_assignment(O)

    E=np.array([0.1,0.9])
    G=np.array([
        [[1.0,0.0,0.0]],
        [[2.0,0.0,0.0]],
    ])
    D=np.zeros((2,2,1,3))
    D[0,1,0,0]=0.3
    D[1,0,0,0]=-0.3

    Et,Gt,Dt=transform_state_properties(E,G,D,result)

    assert np.allclose(Et,[0.9,0.1])
    assert np.allclose(Gt[:,0,0],[2.0,1.0])
    assert Dt[0,1,0,0] == pytest.approx(0.3)
    assert Dt[1,0,0,0] == pytest.approx(-0.3)


def test_ambiguous_assignment_is_detected():
    s=1/np.sqrt(2)
    O=np.array([[s,s],[s,-s]])

    result=maximum_overlap_assignment(
        O,
        minimum_overlap=0.4,
        minimum_score_margin=0.01,
    )

    assert result.ambiguous
    assert any("score margin" in r for r in result.reasons)


def test_energy_clusters_and_subspace_singular_values():
    clusters=energy_degeneracy_clusters(
        np.array([0.0,1e-5,0.1,0.10002]),
        tolerance=5e-5,
    )
    assert clusters == ((0,1),(2,3))

    O=np.eye(3)
    s=subspace_overlap_singular_values(O,(0,1),(0,1))
    assert np.allclose(s,[1.0,1.0])
