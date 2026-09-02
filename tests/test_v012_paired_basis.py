import numpy as np

from gaussian_dynamics.paired_basis_management_v12 import (
    spinor_wavefunction_norm,
    project_spinor_coefficients_to_subset,
    prune_nuclear_gaussian_pairs,
)


def test_paired_projection_preserves_both_electronic_components():
    S=np.array([
        [1.0,0.999999999],
        [0.999999999,1.0],
    ],complex)
    C=np.array([
        [0.8+0.1j,0.2],
        [0.7-0.1j,-0.1j],
    ],complex)

    Cnew,loss=project_spinor_coefficients_to_subset(
        C,S,[0]
    )
    assert Cnew.shape==(1,2)
    assert loss>=0.0
    assert spinor_wavefunction_norm(Cnew,S[np.ix_([0],[0])]) <= spinor_wavefunction_norm(C,S)+1e-12


def test_paired_pruning_removes_whole_redundant_gaussian():
    S=np.array([
        [1.0,0.99999999995],
        [0.99999999995,1.0],
    ],complex)
    C=np.array([
        [0.6,0.2j],
        [0.6,-0.2j],
    ],complex)

    result=prune_nuclear_gaussian_pairs(
        C,S,
        condition_limit=1e8,
        max_projection_loss=1.0,
    )

    assert len(result.keep)==1
    assert result.coefficients_matrix.shape==(1,2)
    assert len(result.removed)==1
