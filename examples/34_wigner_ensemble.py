import numpy as np

from gaussian_dynamics import (
    CIPassageConfig,
    sample_gaussian_wigner,
    gaussian_wigner_covariances,
)

config=CIPassageConfig()
ensemble=sample_gaussian_wigner(
    config.q_array(),
    config.p_array(),
    config.A_matrix(),
    nsamples=10000,
    seed=12345,
)

cov_q,cov_p=gaussian_wigner_covariances(config.A_matrix())

print("Gaussian Wigner initial-condition sampling")
print("------------------------------------------")
print("target q mean:",config.q_array())
print("sample q mean:",np.mean(ensemble.q,axis=0))
print("\ntarget p mean:",config.p_array())
print("sample p mean:",np.mean(ensemble.p,axis=0))
print("\ntarget Cov(q):\n",cov_q)
print("sample Cov(q):\n",np.cov(ensemble.q,rowvar=False))
print("\ntarget Cov(p):\n",cov_p)
print("sample Cov(p):\n",np.cov(ensemble.p,rowvar=False))
