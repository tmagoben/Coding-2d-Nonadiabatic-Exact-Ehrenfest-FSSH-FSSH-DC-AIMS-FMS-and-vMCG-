"""Independent one-dimensional exact-grid reference for v0.22.0 SOC models."""

from dataclasses import asdict, dataclass
import numpy as np

from .matrix_invariants_v213 import hermiticity_residual_v213
from .soc_admission_v221 import (
    require_soc_symmetry_contract_v221,
    soc_symmetry_contract_from_provider_v221,
)


def _positive_dx_v221(dx):
    dx = float(dx)
    if not np.isfinite(dx) or dx <= 0.0:
        raise ValueError("grid spacing must be finite and positive.")
    return dx


def _validate_uniform_grid_v220(x):
    x = np.asarray(x, dtype=float)
    if x.ndim != 1 or len(x) < 8 or not np.all(np.isfinite(x)):
        raise ValueError("spinor grid must contain at least eight finite points.")
    differences = np.diff(x)
    dx = float(differences[0])
    if dx <= 0.0 or not np.allclose(differences, dx, rtol=0.0, atol=1.0e-13):
        raise ValueError("spinor grid must be uniformly increasing.")
    return x, dx


def normalize_spinor_grid_v220(psi, dx):
    psi = np.asarray(psi, dtype=complex)
    dx = _positive_dx_v221(dx)
    if psi.ndim != 2 or not np.all(np.isfinite(psi)):
        raise ValueError("spinor wavefunction must be a finite (nstate,nx) array.")
    norm = float(np.real(np.vdot(psi, psi)) * dx)
    if not np.isfinite(norm) or norm <= 0.0:
        raise ValueError("spinor wavefunction must have finite positive norm.")
    return psi / np.sqrt(norm)


def initial_gaussian_spinor_v220(
    x,
    electronic_vector,
    *,
    center=-1.0,
    momentum=1.2,
    width=0.7,
):
    x, dx = _validate_uniform_grid_v220(x)
    vector = np.asarray(electronic_vector, dtype=complex)
    if vector.ndim != 1 or len(vector) < 1 or not np.all(np.isfinite(vector)):
        raise ValueError("initial electronic vector must be finite and one-dimensional.")
    if np.linalg.norm(vector) == 0.0:
        raise ValueError("initial electronic vector cannot be zero.")
    if not np.isfinite(center) or not np.isfinite(momentum):
        raise ValueError("initial Gaussian center and momentum must be finite.")
    if not np.isfinite(width) or float(width) <= 0.0:
        raise ValueError("initial Gaussian width must be finite and positive.")
    vector = vector / np.linalg.norm(vector)
    nuclear = np.exp(
        -0.5 * float(width) * (x - float(center)) ** 2
        + 1j * float(momentum) * (x - float(center))
    )
    return normalize_spinor_grid_v220(vector[:, None] * nuclear[None, :], dx)


def _fixed_frame_grid_data_v221(provider, x):
    provenance = provider.provenance.validate()
    if not provenance.model_space.fixed_frame:
        raise ValueError(
            "the exact-grid reference requires a fixed electronic frame; "
            "moving-frame providers must first be transported to a global frame."
        )
    contract = soc_symmetry_contract_from_provider_v221(provider)
    require_soc_symmetry_contract_v221(
        provenance.model_space,
        contract,
        provenance=provenance,
    )
    snapshots = tuple(
        provider.evaluate_snapshot(np.asarray([coordinate], dtype=float)).validate()
        for coordinate in x
    )
    potential_values = np.asarray(
        [snapshot.point.H for snapshot in snapshots], dtype=complex
    )
    for matrix in potential_values:
        if hermiticity_residual_v213(matrix) > 1.0e-12:
            raise ValueError("exact-grid potential must be Hermitian at every point.")
    if not np.all(np.isfinite(potential_values)):
        raise ValueError("exact-grid potential contains non-finite data.")
    masses = []
    for snapshot in snapshots:
        mass_matrix = np.asarray(snapshot.point.mass_matrix_q_au, dtype=float)
        if mass_matrix.shape != (1, 1):
            raise ValueError("the one-dimensional grid requires a scalar mass matrix.")
        masses.append(float(mass_matrix[0, 0]))
    if not masses or not np.all(np.isfinite(masses)) or min(masses) <= 0.0:
        raise ValueError("exact-grid mass must be finite and positive.")
    mass = masses[0]
    if max(abs(np.asarray(masses) - mass), default=0.0) > 1.0e-12 * max(
        abs(mass), 1.0
    ):
        raise ValueError("the exact-grid reference requires a constant scalar mass.")
    return potential_values, mass, contract


def _potential_values_v220(provider, x):
    return _fixed_frame_grid_data_v221(provider, x)[0]


def _potential_propagator_v220(potential_values, dt_fraction):
    potential_values = np.asarray(potential_values, dtype=complex)
    if potential_values.ndim != 3 or not np.all(np.isfinite(potential_values)):
        raise ValueError("matrix potential must be a finite (nx,nstate,nstate) array.")
    if potential_values.shape[1] != potential_values.shape[2]:
        raise ValueError("matrix potential must be square in electronic space.")
    for matrix in potential_values:
        if hermiticity_residual_v213(matrix) > 1.0e-12:
            raise ValueError("matrix potential must be Hermitian.")
    energies, vectors = np.linalg.eigh(potential_values)
    phases = np.exp(-1j * float(dt_fraction) * energies)
    return np.einsum(
        "xik,xk,xjk->xij", vectors, phases, vectors.conj(), optimize=True
    )


def _apply_pointwise_v220(propagator, psi):
    return np.einsum("xij,jx->ix", propagator, psi, optimize=True)


def _kinetic_phase_v221(nx, dx, dt, mass):
    k = 2.0 * np.pi * np.fft.fftfreq(int(nx), d=_positive_dx_v221(dx))
    return np.exp(-0.5j * float(dt) * k**2 / float(mass))


def _precomputed_split_step_v221(psi, half_potential, kinetic_phase):
    psi = _apply_pointwise_v220(half_potential, psi)
    psi = np.fft.ifft(
        kinetic_phase[None, :] * np.fft.fft(psi, axis=1), axis=1
    )
    return _apply_pointwise_v220(half_potential, psi)


def spinor_split_operator_step_v220(psi, dx, dt, mass, potential_values):
    """One second-order Strang step for a matrix-valued potential."""
    psi = np.asarray(psi, dtype=complex)
    potential_values = np.asarray(potential_values, dtype=complex)
    dx = _positive_dx_v221(dx)
    if psi.ndim != 2 or potential_values.shape != (psi.shape[1], psi.shape[0], psi.shape[0]):
        raise ValueError("spinor wavefunction and matrix potential dimensions differ.")
    if not np.all(np.isfinite(psi)) or not np.all(np.isfinite(potential_values)):
        raise ValueError("spinor wavefunction and matrix potential must be finite.")
    if not np.isfinite(dt) or float(dt) == 0.0:
        raise ValueError("spinor propagation dt must be finite and nonzero.")
    if not np.isfinite(mass) or float(mass) <= 0.0:
        raise ValueError("spinor propagation mass must be finite and positive.")
    half_potential = _potential_propagator_v220(potential_values, 0.5 * float(dt))
    kinetic_phase = _kinetic_phase_v221(psi.shape[1], dx, dt, mass)
    return _precomputed_split_step_v221(psi, half_potential, kinetic_phase)


def spinor_grid_energy_v220(psi, dx, mass, potential_values):
    psi = np.asarray(psi, dtype=complex)
    potential_values = np.asarray(potential_values, dtype=complex)
    dx = _positive_dx_v221(dx)
    if psi.ndim != 2 or potential_values.shape != (psi.shape[1], psi.shape[0], psi.shape[0]):
        raise ValueError("spinor and potential dimensions differ in the energy audit.")
    if not np.all(np.isfinite(psi)) or not np.all(np.isfinite(potential_values)):
        raise ValueError("spinor and potential must be finite in the energy audit.")
    if not np.isfinite(mass) or float(mass) <= 0.0:
        raise ValueError("energy audit mass must be finite and positive.")
    k = 2.0 * np.pi * np.fft.fftfreq(psi.shape[1], d=dx)
    second = np.fft.ifft(-(k**2)[None, :] * np.fft.fft(psi, axis=1), axis=1)
    kinetic_psi = -second / (2.0 * float(mass))
    potential_psi = np.einsum("xij,jx->ix", potential_values, psi, optimize=True)
    return float(np.real(np.vdot(psi, kinetic_psi + potential_psi)) * dx)


def spinor_grid_projector_population_v220(psi, dx, projector):
    psi = np.asarray(psi, dtype=complex)
    projector = np.asarray(projector, dtype=complex)
    dx = _positive_dx_v221(dx)
    if psi.ndim != 2 or projector.shape != (psi.shape[0], psi.shape[0]):
        raise ValueError("spinor and projector dimensions differ.")
    if not np.all(np.isfinite(psi)) or not np.all(np.isfinite(projector)):
        raise ValueError("spinor and projector must be finite.")
    if max(
        float(np.linalg.norm(projector - projector.conj().T)),
        float(np.linalg.norm(projector @ projector - projector)),
    ) > 1.0e-10:
        raise ValueError("population operator must be a Hermitian projector.")
    density = psi @ psi.conj().T * dx
    trace = float(np.real(np.trace(density)))
    if trace <= 0.0 or not np.isfinite(trace):
        raise ValueError("spinor grid state has invalid norm.")
    return float(np.real(np.trace(density @ projector)) / trace)


def phase_aligned_spinor_grid_error_v220(reference, candidate, dx):
    reference = np.asarray(reference, dtype=complex)
    candidate = np.asarray(candidate, dtype=complex)
    if reference.shape != candidate.shape:
        raise ValueError("spinor grid states have incompatible shapes.")
    dx = _positive_dx_v221(dx)
    if not np.all(np.isfinite(reference)) or not np.all(np.isfinite(candidate)):
        raise ValueError("spinor grid states must be finite.")
    overlap = np.vdot(reference, candidate) * dx
    phase = 1.0 + 0.0j if abs(overlap) < 1.0e-30 else np.exp(-1j * np.angle(overlap))
    squared = float(
        np.real(
            np.vdot(
                phase * candidate - reference,
                phase * candidate - reference,
            )
        )
        * dx
    )
    return float(np.sqrt(max(squared, 0.0)))


@dataclass(frozen=True)
class SpinorGridSettingsV220:
    dt: float = 0.04
    steps: int = 100
    store_every: int = 10

    def validate(self):
        if not np.isfinite(self.dt) or float(self.dt) == 0.0:
            raise ValueError("grid dt must be finite and nonzero.")
        if int(self.steps) != self.steps or int(self.steps) < 0:
            raise ValueError("grid steps must be a nonnegative integer.")
        if int(self.store_every) != self.store_every or int(self.store_every) < 1:
            raise ValueError("grid store_every must be a positive integer.")
        return self


def run_spinor_exact_grid_v220(
    provider,
    x,
    psi0,
    *,
    settings=SpinorGridSettingsV220(),
):
    """Propagate a fixed-frame SOC spinor without Gaussian-dynamics routines."""
    settings = settings.validate()
    x, dx = _validate_uniform_grid_v220(x)
    psi = normalize_spinor_grid_v220(psi0, dx)
    if psi.shape != (provider.provenance.model_space.nstate, len(x)):
        raise ValueError("initial spinor dimension differs from the provider model space.")
    potential_values, mass, symmetry_contract = _fixed_frame_grid_data_v221(
        provider, x
    )
    projectors = {
        name: np.asarray(projector, dtype=complex)
        for name, projector in symmetry_contract.projectors.items()
    }
    half_potential = _potential_propagator_v220(
        potential_values, 0.5 * float(settings.dt)
    )
    kinetic_phase = _kinetic_phase_v221(
        len(x), dx, settings.dt, mass
    )
    times = []
    states = []
    norms = []
    energies = []
    populations = {name: [] for name in projectors}

    def record(step):
        times.append(float(step) * float(settings.dt))
        states.append(psi.copy())
        norms.append(float(np.real(np.vdot(psi, psi)) * dx))
        energies.append(spinor_grid_energy_v220(psi, dx, mass, potential_values))
        for name, projector in projectors.items():
            populations[name].append(
                spinor_grid_projector_population_v220(psi, dx, projector)
            )

    record(0)
    for step in range(1, int(settings.steps) + 1):
        psi = _precomputed_split_step_v221(
            psi, half_potential, kinetic_phase
        )
        if step % int(settings.store_every) == 0 or step == int(settings.steps):
            record(step)
    return {
        "x": x.copy(),
        "dx": dx,
        "time": np.asarray(times),
        "psi": np.asarray(states),
        "norm": np.asarray(norms),
        "energy": np.asarray(energies),
        "populations": {name: np.asarray(values) for name, values in populations.items()},
        "maximum_norm_drift": float(max(abs(np.asarray(norms) - norms[0]))),
        "maximum_energy_drift": float(max(abs(np.asarray(energies) - energies[0]))),
        "settings": asdict(settings),
        "provider_fingerprint": provider.provenance.fingerprint(),
        "release_path": "v0.22.1",
        "fixed_frame_certified": True,
        "constant_mass_certified": True,
    }
