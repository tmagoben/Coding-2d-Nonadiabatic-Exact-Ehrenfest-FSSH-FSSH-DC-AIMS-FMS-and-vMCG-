import numpy as np

from .gaussian import frozen_gaussian


def initial_heller_parameters(q0, p0, sigma):
    """Return q, p, complex A, complex gamma for a normalized Heller packet.

    The convention is
        psi = exp{i[ A(x-q)^2/2 + p(x-q) + gamma ]}

    and the normalized initial Gaussian is
        exp[-(x-q)^2/(4 sigma^2)].
    """
    if sigma <= 0.0:
        raise ValueError("sigma must be positive.")

    A0 = 1j / (2.0 * sigma**2)
    N0 = (1.0 / (2.0 * np.pi * sigma**2)) ** 0.25
    gamma0 = -1j * np.log(N0)

    return float(q0), float(p0), complex(A0), complex(gamma0)


def heller_wavefunction(x, q, p, A, gamma):
    """Reconstruct the Heller wavepacket."""
    y = np.asarray(x, dtype=float) - q
    return np.exp(
        1j * (
            0.5 * A * y**2
            + p * y
            + gamma
        )
    )


def _heller_rhs(state, mass, potential, gradient, hessian):
    q, p, A, gamma = state

    dq = p / mass
    dp = -gradient(q)
    dA = -(A**2) / mass - hessian(q)
    dgamma = p**2 / (2.0 * mass) - potential(q) + 0.5j * A / mass

    return np.array([dq, dp, dA, dgamma], dtype=complex)


def _rk4_complex_step(state, dt, rhs):
    k1 = rhs(state)
    k2 = rhs(state + 0.5 * dt * k1)
    k3 = rhs(state + 0.5 * dt * k2)
    k4 = rhs(state + dt * k3)
    return state + dt * (k1 + 2*k2 + 2*k3 + k4) / 6.0


def run_thawed_gaussian(
    q0,
    p0,
    sigma,
    mass,
    potential,
    gradient,
    hessian,
    dt=0.002,
    steps=1000,
    x=None,
    store_every=1,
):
    """Propagate Heller's thawed Gaussian approximation with RK4."""
    q, p, A, gamma = initial_heller_parameters(q0, p0, sigma)
    state = np.array([q, p, A, gamma], dtype=complex)

    times = []
    parameters = []
    wavefunctions = []

    rhs = lambda s: _heller_rhs(s, mass, potential, gradient, hessian)

    def record(step):
        times.append(step * dt)
        parameters.append(state.copy())
        if x is not None:
            wavefunctions.append(
                heller_wavefunction(x, state[0].real, state[1].real, state[2], state[3])
            )

    record(0)

    for step in range(1, steps + 1):
        state = _rk4_complex_step(state, dt, rhs)

        # q and p are theoretically real for a real potential.
        state[0] = state[0].real
        state[1] = state[1].real

        if step % store_every == 0:
            record(step)

    result = {
        "time": np.asarray(times),
        "parameters": np.asarray(parameters),
        "q": np.asarray(parameters)[:, 0].real,
        "p": np.asarray(parameters)[:, 1].real,
        "A": np.asarray(parameters)[:, 2],
        "gamma": np.asarray(parameters)[:, 3],
    }

    if x is not None:
        result["psi"] = np.asarray(wavefunctions)

    return result


def run_frozen_gaussian(
    q0,
    p0,
    alpha,
    mass,
    potential,
    gradient,
    dt=0.002,
    steps=1000,
    x=None,
    store_every=1,
):
    """Trajectory-guided frozen Gaussian with a fixed width and classical action."""
    q = float(q0)
    p = float(p0)
    action = 0.0

    times = []
    qs = []
    ps = []
    actions = []
    wavefunctions = []

    def record(step):
        times.append(step * dt)
        qs.append(q)
        ps.append(p)
        actions.append(action)
        if x is not None:
            wavefunctions.append(
                frozen_gaussian(x, q, p, alpha) * np.exp(1j * action)
            )

    record(0)

    for step in range(1, steps + 1):
        force = -gradient(q)
        p_half = p + 0.5 * dt * force
        q_new = q + dt * p_half / mass

        # Midpoint action update.
        q_mid = 0.5 * (q + q_new)
        p_mid = p_half
        lagrangian = p_mid**2 / (2.0 * mass) - potential(q_mid)
        action_new = action + dt * lagrangian

        force_new = -gradient(q_new)
        p_new = p_half + 0.5 * dt * force_new

        q, p, action = q_new, p_new, action_new

        if step % store_every == 0:
            record(step)

    result = {
        "time": np.asarray(times),
        "q": np.asarray(qs),
        "p": np.asarray(ps),
        "action": np.asarray(actions),
    }

    if x is not None:
        result["psi"] = np.asarray(wavefunctions)

    return result
