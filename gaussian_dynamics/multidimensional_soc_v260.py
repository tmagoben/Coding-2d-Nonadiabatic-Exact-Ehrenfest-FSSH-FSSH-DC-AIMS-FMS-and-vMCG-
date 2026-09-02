"""Reference-first multidimensional CI+SOC models and exact-grid dynamics.

This module is deliberately independent of the Gaussian TDVP implementation.  It
defines fixed-frame quadratic matrix Hamiltonians and a two-dimensional periodic
FFT/Strang propagator that can therefore act as an external numerical oracle for
the v0.26.0 variational code.

Coordinate convention
---------------------

For ``R in R**D`` the nuclear Hamiltonian is

    H = -1/2 nabla.T @ M^{-1} @ nabla * I + V(R)

with

    V(R) = H0 + sum_a H1[a] R[a]
                  + sum_ab H2[a,b] R[a] R[b].

``H2`` is symmetric in its coordinate indices; no implicit factor of one half is
used.  The electronic frame is global and diabatic.  All matrix coefficients are
Hermitian, and all admitted spin models contain complete declared manifolds.
"""

from dataclasses import asdict, dataclass, field
import hashlib
import json

import numpy as np


MULTIDIMENSIONAL_SOC_MODEL_SCHEMA_V260 = "gnd-multidimensional-soc-model-v0.26.0"
EXACT_GRID_SCHEMA_V260 = "gnd-exact-grid-ci-soc-trajectory-v0.26.0"
QUADRATIC_CONVENTION_V260 = (
    "V(R)=H0+sum_a H1[a]R[a]+sum_ab H2[a,b]R[a]R[b], with H2 coordinate-symmetric"
)
KINETIC_CONVENTION_V260 = "T=-1/2 nabla^T M^{-1} nabla in atomic units"
ELECTRONIC_FRAME_CONVENTION_V260 = "single fixed global diabatic electronic frame"
GRID_INTEGRATOR_V260 = "second-order Strang split operator with unitary FFT kinetic step"


def _canonical_v260(value):
    if isinstance(value, np.generic):
        return _canonical_v260(value.item())
    if isinstance(value, complex):
        return [float(value.real), float(value.imag)]
    if isinstance(value, np.ndarray):
        return _canonical_v260(value.tolist())
    if isinstance(value, dict):
        return {
            str(key): _canonical_v260(item)
            for key, item in sorted(value.items(), key=lambda pair: str(pair[0]))
        }
    if isinstance(value, (tuple, list)):
        return [_canonical_v260(item) for item in value]
    if isinstance(value, float):
        if not np.isfinite(value):
            raise ValueError("v0.26.0 canonical data cannot contain non-finite values.")
        return float(value)
    if value is None or isinstance(value, (str, int, bool)):
        return value
    raise TypeError(f"unsupported v0.26.0 canonical value: {type(value).__name__}")


def _sha256_v260(value):
    payload = json.dumps(
        _canonical_v260(value),
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _complex_pairs_v260(value):
    array = np.asarray(value, dtype=complex)
    return np.stack((array.real, array.imag), axis=-1).tolist()


def _scaled_norm_v260(left, right):
    left = np.asarray(left)
    right = np.asarray(right)
    if left.shape != right.shape:
        return float("inf")
    scale = max(float(np.linalg.norm(left)), float(np.linalg.norm(right)), 1.0)
    return float(np.linalg.norm(left - right) / scale)


def _validate_projectors_v260(projectors, nstate, tolerance):
    result = {}
    for name, projector in dict(projectors).items():
        name = str(name)
        if not name:
            raise ValueError("projector names cannot be empty.")
        projector = np.asarray(projector, dtype=complex)
        if projector.shape != (nstate, nstate) or not np.all(np.isfinite(projector)):
            raise ValueError(f"projector {name!r} has an invalid shape or values.")
        if _scaled_norm_v260(projector, projector.conj().T) > tolerance:
            raise ValueError(f"projector {name!r} is not Hermitian.")
        if _scaled_norm_v260(projector @ projector, projector) > tolerance:
            raise ValueError(f"projector {name!r} is not idempotent.")
        result[name] = projector.copy()
    if result:
        total = sum(result.values(), np.zeros((nstate, nstate), dtype=complex))
        if _scaled_norm_v260(total, np.eye(nstate)) > tolerance:
            raise ValueError("declared electronic projectors must resolve the identity.")
        names = tuple(result)
        for i, left in enumerate(names):
            for right in names[i + 1 :]:
                if float(np.linalg.norm(result[left] @ result[right])) > tolerance:
                    raise ValueError("declared electronic projectors must be orthogonal.")
    return result


@dataclass(frozen=True)
class QuadraticSpinHamiltonianNDV260:
    """Fixed-frame multidimensional Hermitian quadratic spin Hamiltonian."""

    mass_matrix_au: np.ndarray
    H0: np.ndarray
    H1: np.ndarray
    H2: np.ndarray
    label: str = "quadratic multidimensional model"
    model_kind: str = "analytic benchmark"
    projectors: dict = field(default_factory=dict)
    complete_spin_manifold: bool = True
    physical_soc: bool = False
    soc_scale_hartree: float = 0.0
    source: dict = field(default_factory=dict)

    def __post_init__(self):
        object.__setattr__(self, "mass_matrix_au", np.asarray(self.mass_matrix_au, dtype=float).copy())
        object.__setattr__(self, "H0", np.asarray(self.H0, dtype=complex).copy())
        object.__setattr__(self, "H1", np.asarray(self.H1, dtype=complex).copy())
        object.__setattr__(self, "H2", np.asarray(self.H2, dtype=complex).copy())
        object.__setattr__(
            self,
            "projectors",
            {str(key): np.asarray(value, dtype=complex).copy() for key, value in dict(self.projectors).items()},
        )
        object.__setattr__(self, "source", dict(self.source))

    @property
    def ndim(self):
        return int(self.mass_matrix_au.shape[0]) if self.mass_matrix_au.ndim == 2 else 0

    @property
    def nstate(self):
        return int(self.H0.shape[0]) if self.H0.ndim == 2 else 0

    @property
    def inverse_mass_matrix_au(self):
        self.validate()
        return np.linalg.inv(self.mass_matrix_au)

    def validate(self, tolerance=2.0e-11):
        tolerance = float(tolerance)
        if not np.isfinite(tolerance) or tolerance <= 0.0:
            raise ValueError("validation tolerance must be finite and positive.")
        if self.mass_matrix_au.ndim != 2 or self.mass_matrix_au.shape[0] < 1:
            raise ValueError("mass matrix must be nonempty and square.")
        ndim = self.mass_matrix_au.shape[0]
        if self.mass_matrix_au.shape != (ndim, ndim):
            raise ValueError("mass matrix must be square.")
        if not np.all(np.isfinite(self.mass_matrix_au)):
            raise ValueError("mass matrix contains non-finite values.")
        if not np.allclose(self.mass_matrix_au, self.mass_matrix_au.T, atol=tolerance, rtol=0.0):
            raise ValueError("mass matrix must be symmetric.")
        if float(np.min(np.linalg.eigvalsh(self.mass_matrix_au))) <= 0.0:
            raise ValueError("mass matrix must be positive definite.")
        if self.H0.ndim != 2 or self.H0.shape[0] < 1 or self.H0.shape[0] != self.H0.shape[1]:
            raise ValueError("H0 must be a nonempty square matrix.")
        nstate = self.H0.shape[0]
        if self.H1.shape != (ndim, nstate, nstate):
            raise ValueError("H1 must have shape (ndim,nstate,nstate).")
        if self.H2.shape != (ndim, ndim, nstate, nstate):
            raise ValueError("H2 must have shape (ndim,ndim,nstate,nstate).")
        if not all(np.all(np.isfinite(item)) for item in (self.H0, self.H1, self.H2)):
            raise ValueError("Hamiltonian coefficients contain non-finite values.")
        matrices = [self.H0, *self.H1.reshape(-1, nstate, nstate), *self.H2.reshape(-1, nstate, nstate)]
        if any(_scaled_norm_v260(matrix, matrix.conj().T) > tolerance for matrix in matrices):
            raise ValueError("every electronic Hamiltonian coefficient must be Hermitian.")
        if _scaled_norm_v260(self.H2, self.H2.swapaxes(0, 1)) > tolerance:
            raise ValueError("H2 must be symmetric in its coordinate indices.")
        if type(self.complete_spin_manifold) is not bool or type(self.physical_soc) is not bool:
            raise TypeError("spin-manifold and physical-SOC flags must be native Booleans.")
        if not self.complete_spin_manifold:
            raise ValueError("v0.26.0 admits only complete declared spin manifolds.")
        if not np.isfinite(float(self.soc_scale_hartree)) or float(self.soc_scale_hartree) < 0.0:
            raise ValueError("SOC scale must be finite and nonnegative.")
        _validate_projectors_v260(self.projectors, nstate, tolerance)
        if not str(self.label) or not str(self.model_kind):
            raise ValueError("model label and kind cannot be empty.")
        return self

    def hamiltonian(self, coordinates):
        self.validate()
        coordinates = np.asarray(coordinates, dtype=float)
        if coordinates.shape == () or coordinates.shape[-1] != self.ndim:
            raise ValueError("coordinates must have final axis ndim.")
        if not np.all(np.isfinite(coordinates)):
            raise ValueError("coordinates contain non-finite values.")
        value = np.broadcast_to(self.H0, coordinates.shape[:-1] + self.H0.shape).copy()
        value += np.einsum("...a,aij->...ij", coordinates, self.H1, optimize=True)
        value += np.einsum("...a,...b,abij->...ij", coordinates, coordinates, self.H2, optimize=True)
        return value

    def derivative(self, coordinates):
        self.validate()
        coordinates = np.asarray(coordinates, dtype=float)
        if coordinates.shape == () or coordinates.shape[-1] != self.ndim:
            raise ValueError("coordinates must have final axis ndim.")
        if not np.all(np.isfinite(coordinates)):
            raise ValueError("coordinates contain non-finite values.")
        value = np.broadcast_to(self.H1, coordinates.shape[:-1] + self.H1.shape).copy()
        value += 2.0 * np.einsum("...b,abij->...aij", coordinates, self.H2, optimize=True)
        return value

    def gauge_transformed(self, unitary):
        self.validate()
        unitary = np.asarray(unitary, dtype=complex)
        if unitary.shape != (self.nstate, self.nstate):
            raise ValueError("gauge unitary has incompatible shape.")
        if _scaled_norm_v260(unitary.conj().T @ unitary, np.eye(self.nstate)) > 2.0e-11:
            raise ValueError("gauge transformation must be unitary.")

        def rotate(matrix):
            return np.einsum("ia,...ij,jb->...ab", unitary.conj(), matrix, unitary, optimize=True)

        return QuadraticSpinHamiltonianNDV260(
            self.mass_matrix_au,
            rotate(self.H0),
            rotate(self.H1),
            rotate(self.H2),
            label=self.label,
            model_kind=self.model_kind,
            projectors={name: rotate(projector) for name, projector in self.projectors.items()},
            complete_spin_manifold=self.complete_spin_manifold,
            physical_soc=self.physical_soc,
            soc_scale_hartree=self.soc_scale_hartree,
            source={**self.source, "constant_gauge_transformed": True},
        ).validate()

    def coordinate_rotated(self, orthogonal):
        """Return the passive model for ``R_new = O R_old``."""

        self.validate()
        orthogonal = np.asarray(orthogonal, dtype=float)
        if orthogonal.shape != (self.ndim, self.ndim):
            raise ValueError("coordinate rotation has an incompatible shape.")
        if _scaled_norm_v260(orthogonal @ orthogonal.T, np.eye(self.ndim)) > 2.0e-11:
            raise ValueError("coordinate rotation must be orthogonal.")
        return QuadraticSpinHamiltonianNDV260(
            orthogonal @ self.mass_matrix_au @ orthogonal.T,
            self.H0,
            np.einsum("ac,cij->aij", orthogonal, self.H1, optimize=True),
            np.einsum("ac,bd,cdij->abij", orthogonal, orthogonal, self.H2, optimize=True),
            label=self.label + " [coordinate rotated]",
            model_kind=self.model_kind,
            projectors=self.projectors,
            complete_spin_manifold=self.complete_spin_manifold,
            physical_soc=self.physical_soc,
            soc_scale_hartree=self.soc_scale_hartree,
            source={**self.source, "orthogonal_coordinate_rotation": orthogonal.tolist()},
        ).validate()

    def as_dict(self):
        self.validate()
        return {
            "schema": MULTIDIMENSIONAL_SOC_MODEL_SCHEMA_V260,
            "label": str(self.label),
            "model_kind": str(self.model_kind),
            "ndim": self.ndim,
            "nstate": self.nstate,
            "mass_matrix_au": self.mass_matrix_au.tolist(),
            "H0": _complex_pairs_v260(self.H0),
            "H1": _complex_pairs_v260(self.H1),
            "H2": _complex_pairs_v260(self.H2),
            "projectors": {name: _complex_pairs_v260(value) for name, value in self.projectors.items()},
            "complete_spin_manifold": self.complete_spin_manifold,
            "physical_soc": self.physical_soc,
            "soc_scale_hartree": float(self.soc_scale_hartree),
            "quadratic_convention": QUADRATIC_CONVENTION_V260,
            "kinetic_convention": KINETIC_CONVENTION_V260,
            "electronic_frame": ELECTRONIC_FRAME_CONVENTION_V260,
            "source": _canonical_v260(self.source),
        }

    def fingerprint(self):
        return _sha256_v260(self.as_dict())


def _pauli_v260():
    identity = np.eye(2, dtype=complex)
    x = np.asarray([[0.0, 1.0], [1.0, 0.0]], dtype=complex)
    y = np.asarray([[0.0, -1.0j], [1.0j, 0.0]], dtype=complex)
    z = np.asarray([[1.0, 0.0], [0.0, -1.0]], dtype=complex)
    return identity, x, y, z


def two_state_ci_soc_model_v260(
    *,
    mass_au=(1200.0, 1000.0),
    kappa=0.015,
    coupling=0.012,
    frequencies=(0.008, 0.007),
    soc_scale=0.0025,
):
    """Two-state linear-vibronic CI with a constant complex SOC gap."""

    identity, sigma_x, sigma_y, sigma_z = _pauli_v260()
    mass_au = np.asarray(mass_au, dtype=float)
    frequencies = np.asarray(frequencies, dtype=float)
    if mass_au.shape != (2,) or frequencies.shape != (2,):
        raise ValueError("the 2D CI benchmark requires two masses and frequencies.")
    H0 = float(soc_scale) * sigma_y
    H1 = np.stack((float(kappa) * sigma_z, float(coupling) * sigma_x))
    H2 = np.zeros((2, 2, 2, 2), dtype=complex)
    H2[0, 0] = 0.5 * frequencies[0] ** 2 * identity
    H2[1, 1] = 0.5 * frequencies[1] ** 2 * identity
    return QuadraticSpinHamiltonianNDV260(
        np.diag(mass_au),
        H0,
        H1,
        H2,
        label="two-state 2D conical-intersection plus complex SOC",
        model_kind="two-state-ci-soc",
        projectors={"diabatic_0": np.diag([1.0, 0.0]), "diabatic_1": np.diag([0.0, 1.0])},
        complete_spin_manifold=True,
        physical_soc=False,
        soc_scale_hartree=abs(float(soc_scale)),
        source={"analytic": True, "soc_role": "constant sigma_y gap term"},
    ).validate()


def kramers_doublet_ci_soc_model_v260(
    *,
    mass_au=(1200.0, 1000.0),
    kappa=0.015,
    coupling=0.012,
    frequencies=(0.008, 0.007),
    soc_scale=0.0025,
):
    """Two orbital states times a complete Kramers doublet.

    The SOC term ``tau_y tensor sigma_z`` is invariant under spin-1/2 time
    reversal and preserves exact Kramers degeneracy.
    """

    identity, sigma_x, sigma_y, sigma_z = _pauli_v260()
    tau_identity, tau_x, tau_y, tau_z = identity, sigma_x, sigma_y, sigma_z
    electronic_identity = np.eye(4, dtype=complex)
    mass_au = np.asarray(mass_au, dtype=float)
    frequencies = np.asarray(frequencies, dtype=float)
    H0 = float(soc_scale) * np.kron(tau_y, sigma_z)
    H1 = np.stack(
        (
            float(kappa) * np.kron(tau_z, identity),
            float(coupling) * np.kron(tau_x, identity),
        )
    )
    H2 = np.zeros((2, 2, 4, 4), dtype=complex)
    H2[0, 0] = 0.5 * frequencies[0] ** 2 * electronic_identity
    H2[1, 1] = 0.5 * frequencies[1] ** 2 * electronic_identity
    return QuadraticSpinHamiltonianNDV260(
        np.diag(mass_au),
        H0,
        H1,
        H2,
        label="2D CI with a complete Kramers-doublet SOC manifold",
        model_kind="kramers-doublet-ci-soc",
        projectors={
            "orbital_0": np.diag([1.0, 1.0, 0.0, 0.0]),
            "orbital_1": np.diag([0.0, 0.0, 1.0, 1.0]),
        },
        complete_spin_manifold=True,
        physical_soc=False,
        soc_scale_hartree=abs(float(soc_scale)),
        source={"analytic": True, "time_reversal_square": -1, "kramers_complete": True},
    ).validate()


def singlet_triplet_ci_soc_model_v260(
    *,
    mass_au=(1200.0, 1000.0),
    kappa=0.015,
    coupling=0.012,
    frequencies=(0.008, 0.007),
    triplet_offset=0.006,
    soc_scale=0.0015,
):
    """Two CI-forming singlets plus a complete three-component triplet."""

    mass_au = np.asarray(mass_au, dtype=float)
    frequencies = np.asarray(frequencies, dtype=float)
    H0 = np.zeros((5, 5), dtype=complex)
    H0[2:, 2:] = float(triplet_offset) * np.eye(3)
    scale = float(soc_scale) / np.sqrt(2.0)
    coupling_matrix = scale * np.asarray(
        [[1.0, 1.0j, 0.0], [0.0, 1.0, 1.0j]], dtype=complex
    )
    H0[:2, 2:] = coupling_matrix
    H0[2:, :2] = coupling_matrix.conj().T
    H1 = np.zeros((2, 5, 5), dtype=complex)
    H1[0, 0, 0] = float(kappa)
    H1[0, 1, 1] = -float(kappa)
    H1[1, 0, 1] = H1[1, 1, 0] = float(coupling)
    H2 = np.zeros((2, 2, 5, 5), dtype=complex)
    H2[0, 0] = 0.5 * frequencies[0] ** 2 * np.eye(5)
    H2[1, 1] = 0.5 * frequencies[1] ** 2 * np.eye(5)
    return QuadraticSpinHamiltonianNDV260(
        np.diag(mass_au),
        H0,
        H1,
        H2,
        label="2D singlet CI coupled by SOC to a complete triplet",
        model_kind="singlet-triplet-ci-soc",
        projectors={
            "singlet": np.diag([1.0, 1.0, 0.0, 0.0, 0.0]),
            "triplet": np.diag([0.0, 0.0, 1.0, 1.0, 1.0]),
        },
        complete_spin_manifold=True,
        physical_soc=False,
        soc_scale_hartree=abs(float(soc_scale)),
        source={"analytic": True, "singlet_count": 2, "triplet_components": (-1, 0, 1)},
    ).validate()


@dataclass(frozen=True)
class UniformGrid2DV260:
    x: np.ndarray
    y: np.ndarray

    def __post_init__(self):
        object.__setattr__(self, "x", np.asarray(self.x, dtype=float).copy())
        object.__setattr__(self, "y", np.asarray(self.y, dtype=float).copy())

    @classmethod
    def from_bounds(cls, x_bounds, y_bounds, shape):
        nx, ny = (int(item) for item in shape)
        if nx < 8 or ny < 8:
            raise ValueError("each exact-grid dimension requires at least eight points.")
        return cls(
            np.linspace(float(x_bounds[0]), float(x_bounds[1]), nx, endpoint=False),
            np.linspace(float(y_bounds[0]), float(y_bounds[1]), ny, endpoint=False),
        ).validate()

    @property
    def shape(self):
        return (len(self.x), len(self.y))

    @property
    def dx(self):
        return float(self.x[1] - self.x[0])

    @property
    def dy(self):
        return float(self.y[1] - self.y[0])

    @property
    def volume_element(self):
        return self.dx * self.dy

    def mesh(self):
        X, Y = np.meshgrid(self.x, self.y, indexing="ij")
        return np.stack((X, Y), axis=-1)

    def validate(self):
        for name, axis in (("x", self.x), ("y", self.y)):
            if axis.ndim != 1 or len(axis) < 8 or not np.all(np.isfinite(axis)):
                raise ValueError(f"{name} grid must contain at least eight finite points.")
            differences = np.diff(axis)
            if differences[0] <= 0.0 or not np.allclose(differences, differences[0], atol=2.0e-13, rtol=0.0):
                raise ValueError(f"{name} grid must be uniformly increasing.")
        return self

    def as_dict(self):
        self.validate()
        return {
            "shape": list(self.shape),
            "x_min": float(self.x[0]),
            "x_periodic_upper": float(self.x[0] + len(self.x) * self.dx),
            "y_min": float(self.y[0]),
            "y_periodic_upper": float(self.y[0] + len(self.y) * self.dy),
            "dx": self.dx,
            "dy": self.dy,
            "boundary_condition": "periodic FFT box",
        }


@dataclass(frozen=True)
class ExactGridSettingsV260:
    dt_au: float = 0.05
    steps: int = 40
    store_every: int = 10

    def validate(self):
        if not np.isfinite(float(self.dt_au)) or float(self.dt_au) == 0.0:
            raise ValueError("exact-grid dt must be finite and nonzero.")
        for name in ("steps", "store_every"):
            value = getattr(self, name)
            if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)):
                raise TypeError(f"{name} must be an integer.")
        if int(self.steps) < 0 or int(self.store_every) < 1:
            raise ValueError("steps must be nonnegative and store_every positive.")
        return self


def normalize_spinor_grid_v260(psi, grid):
    grid = grid.validate()
    psi = np.asarray(psi, dtype=complex)
    if psi.ndim != 3 or psi.shape[1:] != grid.shape or not np.all(np.isfinite(psi)):
        raise ValueError("spinor grid must have shape (nstate,nx,ny) and finite values.")
    norm = float(np.real(np.vdot(psi, psi)) * grid.volume_element)
    if not np.isfinite(norm) or norm <= 0.0:
        raise ValueError("spinor grid norm must be finite and positive.")
    return psi / np.sqrt(norm)


def initial_gaussian_spinor_2d_v260(
    grid,
    electronic_vector,
    *,
    center=(-2.0, 0.0),
    momentum=(7.0, 0.0),
    widths=(1.0, 1.0),
    chirps=(0.0, 0.0),
):
    """Construct a normalized separable Gaussian spinor without TDVP imports."""

    grid = grid.validate()
    electronic_vector = np.asarray(electronic_vector, dtype=complex)
    center = np.asarray(center, dtype=float)
    momentum = np.asarray(momentum, dtype=float)
    widths = np.asarray(widths, dtype=float)
    chirps = np.asarray(chirps, dtype=float)
    if electronic_vector.ndim != 1 or len(electronic_vector) < 1 or not np.all(np.isfinite(electronic_vector)):
        raise ValueError("electronic vector must be finite and one-dimensional.")
    if float(np.linalg.norm(electronic_vector)) == 0.0:
        raise ValueError("electronic vector cannot be zero.")
    if any(item.shape != (2,) for item in (center, momentum, widths, chirps)):
        raise ValueError("2D Gaussian parameters must contain exactly two values.")
    if not all(np.all(np.isfinite(item)) for item in (center, momentum, widths, chirps)):
        raise ValueError("2D Gaussian parameters must be finite.")
    if np.min(widths) <= 0.0:
        raise ValueError("Gaussian widths must be positive.")
    displacement = grid.mesh() - center
    exponent = np.sum(
        -0.5 * widths * displacement**2
        + 0.5j * chirps * displacement**2
        + 1.0j * momentum * displacement,
        axis=-1,
    )
    nuclear = np.exp(exponent)
    psi = electronic_vector[:, None, None] * nuclear[None, :, :]
    return normalize_spinor_grid_v260(psi, grid)


def _potential_grid_v260(model, grid):
    model = model.validate()
    grid = grid.validate()
    if model.ndim != 2:
        raise ValueError("the v0.26.0 exact-grid oracle is two-dimensional.")
    values = model.hamiltonian(grid.mesh())
    if values.shape != grid.shape + (model.nstate, model.nstate):
        raise AssertionError("model potential grid has an unexpected shape.")
    return values


def _potential_propagator_v260(potential, dt_fraction):
    potential = np.asarray(potential, dtype=complex)
    if potential.ndim != 4 or potential.shape[-1] != potential.shape[-2]:
        raise ValueError("potential must have shape (nx,ny,nstate,nstate).")
    energies, vectors = np.linalg.eigh(potential)
    phases = np.exp(-1.0j * float(dt_fraction) * energies)
    return np.einsum("...ik,...k,...jk->...ij", vectors, phases, vectors.conj(), optimize=True)


def _apply_potential_v260(propagator, psi):
    return np.einsum("xyab,bxy->axy", propagator, psi, optimize=True)


def _kinetic_phase_v260(model, grid, dt):
    inverse_mass = model.inverse_mass_matrix_au
    kx = 2.0 * np.pi * np.fft.fftfreq(len(grid.x), d=grid.dx)
    ky = 2.0 * np.pi * np.fft.fftfreq(len(grid.y), d=grid.dy)
    KX, KY = np.meshgrid(kx, ky, indexing="ij")
    wavevectors = np.stack((KX, KY), axis=-1)
    kinetic_energy = 0.5 * np.einsum("...a,ab,...b->...", wavevectors, inverse_mass, wavevectors, optimize=True)
    return np.exp(-1.0j * float(dt) * kinetic_energy)


def exact_grid_split_step_v260(psi, model, grid, dt, *, half_potential=None, kinetic_phase=None):
    """One reversible unitary 2D matrix-valued Strang step."""

    model = model.validate()
    grid = grid.validate()
    psi = np.asarray(psi, dtype=complex)
    if psi.shape != (model.nstate,) + grid.shape or not np.all(np.isfinite(psi)):
        raise ValueError("spinor shape differs from model/grid dimensions.")
    if not np.isfinite(float(dt)) or float(dt) == 0.0:
        raise ValueError("split-step dt must be finite and nonzero.")
    if half_potential is None:
        half_potential = _potential_propagator_v260(_potential_grid_v260(model, grid), 0.5 * float(dt))
    if kinetic_phase is None:
        kinetic_phase = _kinetic_phase_v260(model, grid, dt)
    psi = _apply_potential_v260(half_potential, psi)
    transformed = np.fft.fftn(psi, axes=(1, 2))
    psi = np.fft.ifftn(kinetic_phase[None, :, :] * transformed, axes=(1, 2))
    return _apply_potential_v260(half_potential, psi)


def exact_grid_norm_v260(psi, grid):
    grid = grid.validate()
    psi = np.asarray(psi, dtype=complex)
    return float(np.real(np.vdot(psi, psi)) * grid.volume_element)


def exact_grid_overlap_v260(bra, ket, grid):
    grid = grid.validate()
    bra = np.asarray(bra, dtype=complex)
    ket = np.asarray(ket, dtype=complex)
    if bra.shape != ket.shape or bra.ndim != 3 or bra.shape[1:] != grid.shape:
        raise ValueError("exact-grid overlap states have incompatible shapes.")
    return complex(np.vdot(bra, ket) * grid.volume_element)


def exact_grid_boundary_probability_v260(psi, grid, *, edge_points=4):
    """Probability in the union of edge strips used to audit periodic wraparound."""

    grid = grid.validate()
    psi = np.asarray(psi, dtype=complex)
    if psi.ndim != 3 or psi.shape[1:] != grid.shape:
        raise ValueError("spinor grid shape is invalid for the boundary audit.")
    edge_points = int(edge_points)
    if edge_points < 1 or 2 * edge_points >= min(grid.shape):
        raise ValueError("edge_points must leave a nonempty grid interior.")
    mask = np.zeros(grid.shape, dtype=bool)
    mask[:edge_points, :] = True
    mask[-edge_points:, :] = True
    mask[:, :edge_points] = True
    mask[:, -edge_points:] = True
    return float(np.sum(np.abs(psi[:, mask]) ** 2) * grid.volume_element)


def exact_grid_reduced_density_v260(psi, grid, *, normalize=True):
    grid = grid.validate()
    psi = np.asarray(psi, dtype=complex)
    if psi.ndim != 3 or psi.shape[1:] != grid.shape:
        raise ValueError("spinor grid shape is invalid for density reduction.")
    density = np.einsum("axy,bxy->ab", psi, psi.conj(), optimize=True) * grid.volume_element
    density = 0.5 * (density + density.conj().T)
    trace = float(np.real(np.trace(density)))
    if normalize:
        if trace <= 0.0 or not np.isfinite(trace):
            raise ValueError("cannot normalize a nonpositive electronic density.")
        density = density / trace
    return density


def exact_grid_energy_v260(psi, model, grid, *, potential=None):
    model = model.validate()
    grid = grid.validate()
    psi = np.asarray(psi, dtype=complex)
    if psi.shape != (model.nstate,) + grid.shape:
        raise ValueError("spinor shape differs from model/grid dimensions.")
    if potential is None:
        potential = _potential_grid_v260(model, grid)
    inverse_mass = model.inverse_mass_matrix_au
    kx = 2.0 * np.pi * np.fft.fftfreq(len(grid.x), d=grid.dx)
    ky = 2.0 * np.pi * np.fft.fftfreq(len(grid.y), d=grid.dy)
    KX, KY = np.meshgrid(kx, ky, indexing="ij")
    wavevectors = np.stack((KX, KY), axis=-1)
    kinetic_energy = 0.5 * np.einsum("...a,ab,...b->...", wavevectors, inverse_mass, wavevectors, optimize=True)
    transformed = np.fft.fftn(psi, axes=(1, 2))
    kinetic_psi = np.fft.ifftn(kinetic_energy[None, :, :] * transformed, axes=(1, 2))
    potential_psi = np.einsum("xyab,bxy->axy", potential, psi, optimize=True)
    value = np.vdot(psi, kinetic_psi + potential_psi) * grid.volume_element
    if abs(float(np.imag(value))) > 2.0e-10:
        raise ValueError("exact-grid energy has a non-negligible imaginary component.")
    return float(np.real(value))


def phase_aligned_grid_error_v260(reference, candidate, grid):
    grid = grid.validate()
    reference = np.asarray(reference, dtype=complex)
    candidate = np.asarray(candidate, dtype=complex)
    if reference.shape != candidate.shape or reference.ndim != 3 or reference.shape[1:] != grid.shape:
        raise ValueError("grid states have incompatible shapes.")
    overlap = np.vdot(reference, candidate) * grid.volume_element
    phase = 1.0 + 0.0j if abs(overlap) < 1.0e-30 else np.exp(-1.0j * np.angle(overlap))
    difference = phase * candidate - reference
    return float(np.sqrt(max(float(np.real(np.vdot(difference, difference)) * grid.volume_element), 0.0)))


@dataclass(frozen=True)
class ExactGridTrajectoryV260:
    model: QuadraticSpinHamiltonianNDV260
    grid: UniformGrid2DV260
    settings: ExactGridSettingsV260
    times_au: np.ndarray
    states: tuple
    norms: np.ndarray
    energies_hartree: np.ndarray
    reduced_densities: np.ndarray
    populations: dict

    def __post_init__(self):
        object.__setattr__(self, "times_au", np.asarray(self.times_au, dtype=float).copy())
        object.__setattr__(self, "states", tuple(np.asarray(item, dtype=complex).copy() for item in self.states))
        object.__setattr__(self, "norms", np.asarray(self.norms, dtype=float).copy())
        object.__setattr__(self, "energies_hartree", np.asarray(self.energies_hartree, dtype=float).copy())
        object.__setattr__(self, "reduced_densities", np.asarray(self.reduced_densities, dtype=complex).copy())
        object.__setattr__(self, "populations", {str(key): np.asarray(value, dtype=float).copy() for key, value in dict(self.populations).items()})

    @property
    def maximum_norm_drift(self):
        return float(np.max(np.abs(self.norms - self.norms[0]))) if len(self.norms) else 0.0

    @property
    def maximum_energy_drift_hartree(self):
        return float(np.max(np.abs(self.energies_hartree - self.energies_hartree[0]))) if len(self.energies_hartree) else 0.0

    @property
    def final_state(self):
        return self.states[-1]

    def validate(self, tolerance=3.0e-10):
        self.model.validate()
        self.grid.validate()
        self.settings.validate()
        count = len(self.times_au)
        if count < 1 or len(self.states) != count:
            raise ValueError("exact-grid trajectory record lengths differ.")
        if self.norms.shape != (count,) or self.energies_hartree.shape != (count,):
            raise ValueError("exact-grid scalar record lengths differ.")
        if self.reduced_densities.shape != (count, self.model.nstate, self.model.nstate):
            raise ValueError("exact-grid density record has an invalid shape.")
        if any(state.shape != (self.model.nstate,) + self.grid.shape for state in self.states):
            raise ValueError("an exact-grid stored state has an invalid shape.")
        arrays = (self.times_au, self.norms, self.energies_hartree, self.reduced_densities, *self.states, *self.populations.values())
        if any(not np.all(np.isfinite(item)) for item in arrays):
            raise ValueError("exact-grid trajectory contains non-finite values.")
        if self.maximum_norm_drift > float(tolerance):
            raise ValueError("exact-grid unitary propagation exceeded its norm gate.")
        for name, values in self.populations.items():
            if values.shape != (count,):
                raise ValueError(f"population record {name!r} has an invalid shape.")
        return self

    def as_dict(self):
        self.validate()
        return {
            "schema": EXACT_GRID_SCHEMA_V260,
            "model_fingerprint": self.model.fingerprint(),
            "grid": self.grid.as_dict(),
            "settings": asdict(self.settings),
            "integrator": GRID_INTEGRATOR_V260,
            "times_au": self.times_au.tolist(),
            "norms": self.norms.tolist(),
            "energies_hartree": self.energies_hartree.tolist(),
            "reduced_densities": _complex_pairs_v260(self.reduced_densities),
            "populations": {name: values.tolist() for name, values in self.populations.items()},
            "maximum_norm_drift": self.maximum_norm_drift,
            "maximum_energy_drift_hartree": self.maximum_energy_drift_hartree,
            "stored_state_sha256": [_sha256_v260(_complex_pairs_v260(state)) for state in self.states],
        }

    def fingerprint(self):
        return _sha256_v260(self.as_dict())


def run_exact_grid_ci_soc_v260(model, grid, psi0, *, settings=ExactGridSettingsV260()):
    """Propagate a 2D SOC spinor with no Gaussian-dynamics dependencies."""

    model = model.validate()
    grid = grid.validate()
    settings = settings.validate()
    if model.ndim != 2:
        raise ValueError("the exact v0.26.0 reference requires a two-dimensional model.")
    psi = normalize_spinor_grid_v260(psi0, grid)
    if psi.shape[0] != model.nstate:
        raise ValueError("initial spinor electronic dimension differs from the model.")
    potential = _potential_grid_v260(model, grid)
    half_potential = _potential_propagator_v260(potential, 0.5 * float(settings.dt_au))
    kinetic_phase = _kinetic_phase_v260(model, grid, settings.dt_au)
    times = []
    states = []
    norms = []
    energies = []
    densities = []
    populations = {name: [] for name in model.projectors}

    def record(step):
        density = exact_grid_reduced_density_v260(psi, grid)
        times.append(float(step) * float(settings.dt_au))
        states.append(psi.copy())
        norms.append(exact_grid_norm_v260(psi, grid))
        energies.append(exact_grid_energy_v260(psi, model, grid, potential=potential))
        densities.append(density)
        for name, projector in model.projectors.items():
            populations[name].append(float(np.real(np.trace(density @ projector))))

    record(0)
    for step in range(1, int(settings.steps) + 1):
        psi = exact_grid_split_step_v260(
            psi,
            model,
            grid,
            settings.dt_au,
            half_potential=half_potential,
            kinetic_phase=kinetic_phase,
        )
        if step % int(settings.store_every) == 0 or step == int(settings.steps):
            record(step)
    return ExactGridTrajectoryV260(
        model=model,
        grid=grid,
        settings=settings,
        times_au=np.asarray(times),
        states=tuple(states),
        norms=np.asarray(norms),
        energies_hartree=np.asarray(energies),
        reduced_densities=np.asarray(densities),
        populations={name: np.asarray(values) for name, values in populations.items()},
    ).validate()


V260_EXACT_GRID_CLAIMS = {
    "two_dimensional_ci_soc_exact_grid_validated": True,
    "complete_kramers_doublet_model_validated": True,
    "complete_singlet_triplet_model_validated": True,
    "complex_soc_and_zero_soc_toggle_validated": True,
    "full_positive_definite_mass_matrix_supported": True,
    "unitary_reversible_split_operator_validated": True,
    "independent_of_gaussian_tdvp_implementation": True,
    "absorbing_boundary_conditions_validated": False,
    "multidimensional_grid_beyond_two_dimensions_validated": False,
    "real_ab_initio_potential_grid_validated": False,
}
