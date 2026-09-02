import numpy as np

from gaussian_dynamics.benchmark_provider_nd import LVC2DGeneralizedProvider
from gaussian_dynamics.local_gaussian_nd import LocalAdiabaticTBF
from gaussian_dynamics.direct_dynamics_nd import run_backend_spawned_gaussians

provider=LVC2DGeneralizedProvider(nuclear_mass_au=20.0)

parent=LocalAdiabaticTBF(
    state=1,
    q=np.array([0.55,0.45]),
    p=np.array([0.6,0.8]),
    A=1.2*np.eye(2),
)

out=run_backend_spawned_gaussians(
    initial_basis=[parent],
    C0=[1.0+0.0j],
    provider=provider,
    dt=0.0002,
    steps=50,
    spawn_threshold=1e-6,
    overlap_block=0.9,
    max_basis=2,
    store_every=5,
)

print("v0.5 gridless backend-driven Gaussian dynamics")
print("spawn events:",out["events"])
print("basis size:",out["basis_size"])
print("final populations:",out["state_populations"][-1])
print("max norm drift:",np.max(np.abs(out["norm"]-1.0)))
