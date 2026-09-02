import numpy as np

from gaussian_dynamics.initial_conditions import (
    gaussian_wigner_covariances,
    sample_gaussian_wigner,
)


def test_wigner_covariance_formula():
    A=np.array([[1.2,0.2],[0.2,0.9]])
    qcov,pcov=gaussian_wigner_covariances(A)

    assert np.allclose(qcov,0.5*np.linalg.inv(A))
    assert np.allclose(pcov,0.5*A)


def test_wigner_sampling_is_seed_reproducible():
    A=np.eye(2)
    a=sample_gaussian_wigner([0,0],[1,-1],A,5,seed=77)
    b=sample_gaussian_wigner([0,0],[1,-1],A,5,seed=77)

    assert np.allclose(a.q,b.q)
    assert np.allclose(a.p,b.p)


def test_wigner_empirical_moments_are_reasonable():
    A=np.array([[1.1,0.1],[0.1,0.8]])
    q0=np.array([0.3,-0.4])
    p0=np.array([0.7,0.2])

    ensemble=sample_gaussian_wigner(q0,p0,A,20000,seed=11)
    qcov,pcov=gaussian_wigner_covariances(A)

    assert np.allclose(np.mean(ensemble.q,axis=0),q0,atol=0.02)
    assert np.allclose(np.mean(ensemble.p,axis=0),p0,atol=0.02)
    assert np.allclose(np.cov(ensemble.q,rowvar=False),qcov,atol=0.02)
    assert np.allclose(np.cov(ensemble.p,rowvar=False),pcov,atol=0.02)
