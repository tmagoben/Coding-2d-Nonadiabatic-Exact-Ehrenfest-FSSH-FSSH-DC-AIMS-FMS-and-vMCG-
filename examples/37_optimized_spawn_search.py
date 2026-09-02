import numpy as np

from gaussian_dynamics.dynamic_graph_aims import DynamicGraphTBF
from gaussian_dynamics.dynamic_gauge_graph import AnalyticCI2DFrameProvider
from gaussian_dynamics.optimized_spawning import generate_spawn_candidates

provider=AnalyticCI2DFrameProvider(nuclear_mass_au=20.0)

parent=DynamicGraphTBF(
    uid=0,
    state=1,
    q=np.array([0.55,0.45]),
    p=np.array([0.6,0.8]),
    A=1.2*np.eye(2),
    node=("seed",0),
)

candidates=generate_spawn_candidates(
    parent,
    target=0,
    provider=provider,
    basis=[parent],
    position_shifts=(0.0,0.05,-0.05),
    width_scales=(0.65,1.0,1.55),
    momentum_directions=("nac","momentum"),
    overlap_block=0.9999,
)

print("Top optimized-spawn-inspired local candidates")
print("---------------------------------------------")
for rank,c in enumerate(candidates[:8],start=1):
    print(
        f"{rank:2d} score={c.score:.6e} "
        f"coupling={c.coupling_proxy:.6e} "
        f"shift={c.position_shift:+.3f}({c.position_direction}) "
        f"width_scale={c.width_scale:.3f} "
        f"momentum={c.momentum_direction} "
        f"|S_nuc|={c.nuclear_overlap:.6f} "
        f"dE={c.energy_residual:.3e}"
    )
