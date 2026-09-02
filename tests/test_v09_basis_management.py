import numpy as np

from gaussian_dynamics.basis_management import (
    overlap_conditioning,
    project_coefficients_to_subset,
    prune_redundant_basis,
    canonical_orthogonalizer,
)


def test_projection_removes_duplicate_basis_with_zero_loss():
    S=np.array([[1.0,1.0],[1.0,1.0]],complex)
    C=np.array([0.5,0.5],complex)
    Cnew,loss=project_coefficients_to_subset(C,S,[0])
    assert abs(loss) < 1e-12
    assert np.allclose(Cnew,[1.0])


def test_pruning_improves_conditioning_with_small_loss():
    eps=1e-10
    S=np.array([[1.0,1.0-eps],[1.0-eps,1.0]],complex)
    C=np.array([0.5,0.5],complex)
    result=prune_redundant_basis(
        C,S,condition_limit=1e8,eigenvalue_floor=1e-8,max_projection_loss=1e-6
    )
    assert len(result.removed)==1
    assert result.condition_after < result.condition_before
    assert result.projection_loss < 1e-6


def test_canonical_orthogonalizer_has_identity_metric():
    S=np.array([[1.0,0.3],[0.3,1.0]],complex)
    X,eig,mask=canonical_orthogonalizer(S)
    assert np.allclose(X.conj().T@S@X,np.eye(2),atol=1e-12)
