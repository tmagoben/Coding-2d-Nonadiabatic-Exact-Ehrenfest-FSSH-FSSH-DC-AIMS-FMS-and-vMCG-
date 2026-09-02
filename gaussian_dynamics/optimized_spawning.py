from dataclasses import dataclass
import numpy as np

from .gaussian_general import (
    gaussian_overlap_general,
    gaussian_cross_centroid,
    real_overlap_saddle_point,
    width_scaled,
)


@dataclass(frozen=True)
class SpawnCandidate:
    target_state: int
    q: np.ndarray
    p: np.ndarray
    A: np.ndarray
    position_direction: str
    position_shift: float
    momentum_direction: str
    width_scale: float
    energy_residual: float
    nuclear_overlap: float
    max_existing_overlap: float
    coupling_proxy: float
    novelty: float
    score: float


def classical_energy(q,p,state,provider):
    point=provider.evaluate(np.asarray(q,float))
    Minv=np.linalg.inv(point.mass_matrix)
    return float(point.energies[state] + 0.5*np.asarray(p,float)@Minv@np.asarray(p,float))


def _unit(v, floor=1e-14):
    v=np.asarray(v,float)
    n=np.linalg.norm(v)
    if n<floor:
        return None
    return v/n


def _candidate_position_directions(parent,target,provider):
    point=provider.evaluate(parent.q)

    directions=[]

    nac=_unit(point.nac[parent.state,target])
    if nac is not None:
        directions.append(("nac",nac))

    force_target=_unit(-point.gradients[target])
    if force_target is not None:
        directions.append(("target_force",force_target))

    momentum=_unit(parent.p)
    if momentum is not None:
        directions.append(("momentum",momentum))

    # Remove numerically duplicate directions up to sign: signs are supplied
    # independently through positive/negative shift values.
    unique=[]
    for name,v in directions:
        if any(abs(np.dot(v,u))>1-1e-10 for _,u in unique):
            continue
        unique.append((name,v))

    return unique


def _momentum_adjustment_directions(q,parent,target,provider,names):
    point=provider.evaluate(q)
    out=[]

    for name in names:
        if name=="nac":
            v=_unit(point.nac[parent.state,target])
        elif name=="momentum":
            v=_unit(parent.p)
        elif name=="target_force":
            v=_unit(-point.gradients[target])
        else:
            raise ValueError(f"unknown momentum adjustment direction {name!r}")

        if v is not None:
            out.append((name,v))

    return out


def energy_conserving_momentum_at_position(
    parent,
    target,
    q_child,
    direction,
    provider,
):
    """Adjust parent momentum along `direction` to match parent total energy."""
    q_child=np.asarray(q_child,float)
    n=_unit(direction)
    if n is None:
        return None

    parent_energy=classical_energy(
        parent.q,parent.p,parent.state,provider
    )

    point=provider.evaluate(q_child)
    Minv=np.linalg.inv(point.mass_matrix)
    p0=np.asarray(parent.p,float)

    # 1/2 (p0+lambda n)^T B (p0+lambda n) + E_b(q_child) = E_parent
    a=float(n@Minv@n)
    b=float(p0@Minv@n)
    c=float(p0@Minv@p0 + 2.0*(point.energies[target]-parent_energy))

    disc=b*b-a*c
    if disc < -1e-13:
        return None
    disc=max(0.0,disc)
    root=np.sqrt(disc)

    lambdas=[
        (-b+root)/a,
        (-b-root)/a,
    ]
    lam=min(lambdas,key=abs)

    p_child=p0+lam*n
    return p_child


def local_spa1_coupling_proxy(parent,target,q_child,p_child,A_child,provider):
    """First-order local derivative-Hamiltonian parent-child coupling proxy.

    At a common adiabatic saddle frame the zeroth-order electronic Hamiltonian is
    diagonal.  The leading off-diagonal Taylor term is

        S_nuc * F_ab(q_s) . (mu-q_s),

    with F_ab = (E_b-E_a) d_ab.

    The magnitude is used only to rank spawn candidates.
    """
    qs=real_overlap_saddle_point(
        parent.q,parent.A,q_child,A_child
    )
    mu=gaussian_cross_centroid(
        parent.q,parent.p,parent.A,
        q_child,p_child,A_child,
    )
    S=gaussian_overlap_general(
        parent.q,parent.p,parent.A,
        q_child,p_child,A_child,
    )

    point=provider.evaluate(qs)
    gap=float(point.energies[target]-point.energies[parent.state])
    F=gap*np.asarray(point.nac[parent.state,target],float)

    coupling=S*(F @ (mu-qs))
    return float(abs(coupling)),float(abs(S))


def maximum_target_overlap(candidate_q,candidate_p,candidate_A,target,basis):
    overlaps=[]
    for existing in basis:
        if existing.state != target:
            continue
        overlaps.append(abs(
            gaussian_overlap_general(
                existing.q,existing.p,existing.A,
                candidate_q,candidate_p,candidate_A,
            )
        ))
    return float(max(overlaps,default=0.0))


def generate_spawn_candidates(
    parent,
    target,
    provider,
    basis,
    position_shifts=(0.0,0.05,0.10,-0.05,-0.10),
    width_scales=(0.65,1.0,1.55),
    momentum_directions=("nac","momentum"),
    overlap_block=0.9995,
    novelty_power=0.5,
    energy_tolerance=1e-10,
):
    """Generate energy-conserving, nonredundant local spawn candidates."""
    position_dirs=_candidate_position_directions(parent,target,provider)

    # Zero displacement is direction independent and should be evaluated once.
    q_specs=[("none",np.zeros_like(parent.q),0.0)]
    for name,direction in position_dirs:
        for shift in position_shifts:
            shift=float(shift)
            if abs(shift)<1e-15:
                continue
            q_specs.append((name,direction,shift))

    candidates=[]

    for pos_name,pos_direction,shift in q_specs:
        q_child=parent.q + shift*pos_direction

        for mom_name,mom_direction in _momentum_adjustment_directions(
            q_child,parent,target,provider,momentum_directions
        ):
            p_child=energy_conserving_momentum_at_position(
                parent,target,q_child,mom_direction,provider
            )
            if p_child is None:
                continue

            for scale in width_scales:
                A_child=width_scaled(parent.A,scale)

                max_overlap=maximum_target_overlap(
                    q_child,p_child,A_child,target,basis
                )
                if max_overlap >= overlap_block:
                    continue

                coupling,nuclear_overlap=local_spa1_coupling_proxy(
                    parent,target,q_child,p_child,A_child,provider
                )

                novelty=max(
                    0.0,
                    1.0-max_overlap**2,
                )**float(novelty_power)

                score=coupling*novelty

                e_parent=classical_energy(
                    parent.q,parent.p,parent.state,provider
                )
                e_child=classical_energy(
                    q_child,p_child,target,provider
                )
                residual=e_child-e_parent

                if abs(residual)>energy_tolerance:
                    continue

                candidates.append(SpawnCandidate(
                    target_state=int(target),
                    q=np.asarray(q_child,float),
                    p=np.asarray(p_child,float),
                    A=np.asarray(A_child,float),
                    position_direction=pos_name,
                    position_shift=float(shift),
                    momentum_direction=mom_name,
                    width_scale=float(scale),
                    energy_residual=float(residual),
                    nuclear_overlap=float(nuclear_overlap),
                    max_existing_overlap=float(max_overlap),
                    coupling_proxy=float(coupling),
                    novelty=float(novelty),
                    score=float(score),
                ))

    candidates.sort(key=lambda c:c.score,reverse=True)
    return candidates


def select_spawn_children(
    parent,
    target,
    provider,
    basis,
    children_per_event=1,
    child_overlap_block=0.985,
    **candidate_kwargs,
):
    """Select top ranked candidates while preventing sibling redundancy."""
    ranked=generate_spawn_candidates(
        parent,target,provider,basis,**candidate_kwargs
    )

    selected=[]
    for candidate in ranked:
        if candidate.score<=0.0:
            continue

        redundant=False
        for other in selected:
            ov=abs(gaussian_overlap_general(
                other.q,other.p,other.A,
                candidate.q,candidate.p,candidate.A,
            ))
            if ov>=child_overlap_block:
                redundant=True
                break

        if redundant:
            continue

        selected.append(candidate)
        if len(selected)>=int(children_per_event):
            break

    return selected
