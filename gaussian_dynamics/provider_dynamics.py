from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class ProviderTBF:
    state: int
    q: float
    p: float
    alpha: float = 1.0


def velocity_verlet_tbf(tbf, provider, mass, dt):
    """One trajectory-guidance step using provider adiabatic gradients."""
    p0=provider.evaluate(tbf.q)
    force=-p0.gradients_q[tbf.state]
    p_half=tbf.p+0.5*dt*force
    q_new=tbf.q+dt*p_half/mass

    p1=provider.evaluate(q_new)
    force_new=-p1.gradients_q[tbf.state]
    p_new=p_half+0.5*dt*force_new

    return ProviderTBF(tbf.state,q_new,p_new,tbf.alpha)


def coupling_indicator(tbf, provider, mass, target=None):
    point=provider.evaluate(tbf.q)
    if target is None:
        if len(point.energies) != 2:
            raise ValueError("target is required for more than two states.")
        target=1-tbf.state
    return abs((tbf.p/mass)*point.nac_q[tbf.state,target])


def energy_conserving_child(tbf, provider, mass, target=None):
    point=provider.evaluate(tbf.q)
    if target is None:
        if len(point.energies) != 2:
            raise ValueError("target is required for more than two states.")
        target=1-tbf.state

    rad=tbf.p**2+2*mass*(point.energies[tbf.state]-point.energies[target])
    if rad < 0:
        return None

    sign=1.0 if tbf.p >= 0 else -1.0
    return ProviderTBF(target,tbf.q,sign*np.sqrt(rad),tbf.alpha)


def provider_grid(provider, q_grid):
    """Sample a provider onto a 1D grid for Gaussian/grid benchmark modules."""
    pts=[provider.evaluate(float(q)) for q in q_grid]
    return {
        "q":np.asarray(q_grid,float),
        "energies":np.asarray([p.energies for p in pts]),
        "gradients_q":np.asarray([p.gradients_q for p in pts]),
        "nac_q":np.asarray([p.nac_q for p in pts]),
        "metadata":[p.metadata for p in pts],
    }
