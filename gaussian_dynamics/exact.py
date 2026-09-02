import numpy as np

from .grids import normalize


def split_operator_step(psi, x, dx, dt, mass, potential_values):
    """Second-order Strang split-operator step."""
    psi = np.asarray(psi, dtype=complex)
    V = np.asarray(potential_values, dtype=float)

    half_V = np.exp(-0.5j * dt * V)
    psi = half_V * psi

    k = 2.0 * np.pi * np.fft.fftfreq(len(x), d=dx)
    kinetic_phase = np.exp(-0.5j * dt * k**2 / mass)

    psi_k = np.fft.fft(psi)
    psi = np.fft.ifft(kinetic_phase * psi_k)

    psi = half_V * psi
    return psi


def energy_expectation(psi, x, dx, mass, potential_values):
    """Grid expectation value of T+V."""
    psi = np.asarray(psi, dtype=complex)
    V = np.asarray(potential_values, dtype=float)

    k = 2.0 * np.pi * np.fft.fftfreq(len(x), d=dx)
    second = np.fft.ifft(-(k**2) * np.fft.fft(psi))
    Tpsi = -second / (2.0 * mass)
    Hpsi = Tpsi + V * psi

    return (np.vdot(psi, Hpsi) * dx).real


def run_split_operator(
    psi0,
    x,
    dx,
    potential,
    mass=1.0,
    dt=0.002,
    steps=1000,
    renormalize=False,
    store_every=1,
):
    """Propagate a 1D wavepacket on a time-independent potential."""
    psi = normalize(np.asarray(psi0, dtype=complex), dx)
    V = np.asarray(potential(x), dtype=float)

    times = []
    states = []
    norms = []
    energies = []

    def record(step):
        times.append(step * dt)
        states.append(psi.copy())
        norms.append((np.vdot(psi, psi) * dx).real)
        energies.append(energy_expectation(psi, x, dx, mass, V))

    record(0)

    for step in range(1, steps + 1):
        psi = split_operator_step(psi, x, dx, dt, mass, V)

        if renormalize:
            psi = normalize(psi, dx)

        if step % store_every == 0:
            record(step)

    return {
        "time": np.asarray(times),
        "psi": np.asarray(states),
        "norm": np.asarray(norms),
        "energy": np.asarray(energies),
    }
