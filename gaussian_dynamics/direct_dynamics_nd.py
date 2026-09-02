import numpy as np

from .gaussian_nd import analytic_overlap_equal_width
from .local_gaussian_nd import (
    LocalAdiabaticTBF,
    local_matrices,
    tbf_guidance,
)


def coefficient_rhs(C, S, H, T):
    return np.linalg.solve(S, -1j*(H @ C) - T @ C)


def nac_indicator(tbf, target, provider):
    point = provider.evaluate(tbf.q)
    velocity = np.linalg.solve(point.mass_matrix_q_au, tbf.p)
    return float(abs(velocity @ point.nac_q[tbf.state, target]))


def energy_conserving_child(tbf, target, provider):
    point = provider.evaluate(tbf.q)
    d = np.asarray(point.nac_q[tbf.state, target], dtype=float)
    dn = np.linalg.norm(d)
    if dn < 1e-14:
        return None

    n = d/dn
    B = np.linalg.inv(point.mass_matrix_q_au)

    delta_E = point.energies[target] - point.energies[tbf.state]

    Acoef = float(n @ B @ n)
    Bcoef = float(tbf.p @ B @ n)

    disc = Bcoef**2 - 2.0*Acoef*delta_E
    if disc < 0.0:
        return None

    root = np.sqrt(disc)
    candidates = [
        (-Bcoef + root)/Acoef,
        (-Bcoef - root)/Acoef,
    ]
    lam = min(candidates, key=abs)

    p_child = tbf.p + lam*n

    return LocalAdiabaticTBF(
        state=int(target),
        q=tbf.q.copy(),
        p=p_child,
        A=tbf.A.copy(),
    )


def phase_space_overlap_magnitude(a, b):
    if a.state != b.state:
        return 0.0
    if not np.allclose(a.A, b.A, atol=1e-12):
        return 0.0
    return float(abs(
        analytic_overlap_equal_width(a.q, a.p, b.q, b.p, a.A)
    ))


def maybe_spawn_once(
    basis,
    provider,
    threshold=1e-4,
    overlap_block=0.85,
):
    """Deterministic first-eligible parent/target spawn."""
    for parent_index, parent in enumerate(basis):
        point = provider.evaluate(parent.q)
        ns = len(point.energies)

        for target in range(ns):
            if target == parent.state:
                continue

            if nac_indicator(parent, target, provider) <= threshold:
                continue

            child = energy_conserving_child(parent, target, provider)
            if child is None:
                continue

            redundant = any(
                existing.state == child.state
                and phase_space_overlap_magnitude(existing, child) >= overlap_block
                for existing in basis
            )
            if redundant:
                continue

            return parent_index, child

    return None, None


def _temporary_basis(states, widths, q, p):
    return [
        LocalAdiabaticTBF(int(states[i]), q[i].copy(), p[i].copy(), widths[i].copy())
        for i in range(len(states))
    ]


def run_backend_spawned_gaussians(
    initial_basis,
    C0,
    provider,
    dt=0.0005,
    steps=100,
    spawn_threshold=1e-4,
    overlap_block=0.85,
    max_basis=8,
    store_every=5,
):
    """Gridless backend-driven spawned Gaussian dynamics.

    Electronic quantities are evaluated at TBF centers and pair centroids through the
    provider. Gaussian matrix elements use the local constant-electronic-quantity
    approximation documented in V05_THEORY.md.
    """
    basis = [b.copy() for b in initial_basis]
    C = np.asarray(C0, dtype=complex).copy()

    if len(C) != len(basis):
        raise ValueError("C0 must contain one coefficient per initial TBF.")

    S, _, _ = local_matrices(basis, provider)
    norm0 = float(np.real(np.vdot(C, S @ C)))
    if norm0 <= 0.0:
        raise ValueError("Initial wavefunction has non-positive norm.")
    C /= np.sqrt(norm0)

    times=[]; norms=[]; sizes=[]; conds=[]; state_pops=[]; events=[]

    def diagnostics(step):
        S, H, T = local_matrices(basis, provider)
        norm = float(np.real(np.vdot(C, S @ C)))

        ns = max(b.state for b in basis) + 1
        pops = np.zeros(ns, dtype=float)
        for state in range(ns):
            idx = [i for i,b in enumerate(basis) if b.state == state]
            if idx:
                block = S[np.ix_(idx,idx)]
                cc = C[idx]
                pops[state] = float(np.real(np.vdot(cc, block @ cc)))
        if norm > 0:
            pops /= norm

        times.append(step*dt)
        norms.append(norm)
        sizes.append(len(basis))
        conds.append(np.linalg.cond(S))
        state_pops.append(pops)

    diagnostics(0)

    for step in range(1, steps+1):
        n = len(basis)
        states = np.array([b.state for b in basis], dtype=int)
        widths = [b.A.copy() for b in basis]
        q = np.array([b.q for b in basis], dtype=float)
        p = np.array([b.p for b in basis], dtype=float)

        def rhs(Cx, qx, px):
            temp = _temporary_basis(states, widths, qx, px)
            S, H, T = local_matrices(temp, provider)
            Cdot = coefficient_rhs(Cx, S, H, T)

            qdot = np.zeros_like(qx)
            pdot = np.zeros_like(px)
            for i, b in enumerate(temp):
                qdot[i], pdot[i] = tbf_guidance(b, provider)
            return Cdot, qdot, pdot

        k1 = rhs(C, q, p)
        k2 = rhs(
            C + 0.5*dt*k1[0],
            q + 0.5*dt*k1[1],
            p + 0.5*dt*k1[2],
        )
        k3 = rhs(
            C + 0.5*dt*k2[0],
            q + 0.5*dt*k2[1],
            p + 0.5*dt*k2[2],
        )
        k4 = rhs(
            C + dt*k3[0],
            q + dt*k3[1],
            p + dt*k3[2],
        )

        C = C + dt*(k1[0] + 2*k2[0] + 2*k3[0] + k4[0])/6.0
        q = q + dt*(k1[1] + 2*k2[1] + 2*k3[1] + k4[1])/6.0
        p = p + dt*(k1[2] + 2*k2[2] + 2*k3[2] + k4[2])/6.0

        for i, b in enumerate(basis):
            b.q = q[i].copy()
            b.p = p[i].copy()

        if len(basis) < max_basis:
            parent_idx, child = maybe_spawn_once(
                basis,
                provider,
                threshold=spawn_threshold,
                overlap_block=overlap_block,
            )

            if child is not None:
                basis.append(child)
                C = np.concatenate([C, [0.0+0.0j]])
                events.append({
                    "step": step,
                    "time": step*dt,
                    "parent_index": int(parent_idx),
                    "new_index": len(basis)-1,
                    "target_state": int(child.state),
                })

        if step % store_every == 0:
            diagnostics(step)

    max_states = max(len(p) for p in state_pops)
    padded = np.zeros((len(state_pops), max_states))
    for i,pop in enumerate(state_pops):
        padded[i,:len(pop)] = pop

    return {
        "time": np.asarray(times),
        "norm": np.asarray(norms),
        "basis_size": np.asarray(sizes),
        "condition_number": np.asarray(conds),
        "state_populations": padded,
        "events": events,
        "final_basis": basis,
        "final_coefficients": C,
    }
