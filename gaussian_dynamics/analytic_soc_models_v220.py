"""Reproducible physical analytic SOC models for v0.22.0.

The singlet--triplet and doublet models deliberately occupy separate electron-number
sectors.  Both use one nuclear coordinate, a fixed spin-diabatic frame, complete
multiplets, analytic physical operator derivatives, and explicit time-reversal
conventions.  They are small enough for independent exact-grid propagation.
"""

from dataclasses import asdict, dataclass
import numpy as np

from .electronic_contract_v213 import (
    ElectronicModelSpaceV213,
    ElectronicOperatorProvenanceV213,
    ElectronicStateDescriptorV213,
    compose_electronic_operator_v213,
)
from .electronic_operator_v21 import ElectronicOperatorSnapshotV21
from .matrix_invariants_v213 import hermiticity_residual_v213
from .soc_admission_v221 import SOCSymmetryContractV221


def _finite_real(name, value, *, positive=False):
    value = float(value)
    if not np.isfinite(value) or (positive and value <= 0.0):
        qualifier = "finite and positive" if positive else "finite"
        raise ValueError(f"{name} must be {qualifier}.")
    return value


def singlet_triplet_time_reversal_matrix_v220():
    """Unitary part of even-electron time reversal in (S,T-1,T0,T+1)."""
    return np.asarray(
        [
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
            [0.0, 0.0, -1.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
        ],
        dtype=complex,
    )


def kramers_time_reversal_matrix_v220():
    """Unitary part of odd-electron time reversal for two ordered doublets."""
    j2 = np.asarray([[0.0, 1.0], [-1.0, 0.0]], dtype=complex)
    return np.kron(np.eye(2, dtype=complex), j2)


def singlet_triplet_projectors_v220():
    singlet = np.diag([1.0, 0.0, 0.0, 0.0]).astype(complex)
    triplet = np.eye(4, dtype=complex) - singlet
    return {"singlet": singlet, "triplet": triplet}


def doublet_root_projectors_v220():
    root_1 = np.diag([1.0, 1.0, 0.0, 0.0]).astype(complex)
    root_2 = np.eye(4, dtype=complex) - root_1
    return {"doublet_1": root_1, "doublet_2": root_2}


@dataclass(frozen=True)
class SOCOperatorComponentsV220:
    q: np.ndarray
    H_spin_free: np.ndarray
    K_spin_free: np.ndarray
    H_soc: np.ndarray
    K_soc: np.ndarray

    def __post_init__(self):
        object.__setattr__(self, "q", np.asarray(self.q, dtype=float).copy())
        object.__setattr__(
            self, "H_spin_free", np.asarray(self.H_spin_free, dtype=complex).copy()
        )
        object.__setattr__(
            self, "K_spin_free", np.asarray(self.K_spin_free, dtype=complex).copy()
        )
        object.__setattr__(self, "H_soc", np.asarray(self.H_soc, dtype=complex).copy())
        object.__setattr__(self, "K_soc", np.asarray(self.K_soc, dtype=complex).copy())

    @property
    def H(self):
        return self.H_spin_free + self.H_soc

    @property
    def K(self):
        return self.K_spin_free + self.K_soc

    def validate(self, tolerance=1.0e-12):
        q = np.asarray(self.q, dtype=float)
        H0 = np.asarray(self.H_spin_free, dtype=complex)
        K0 = np.asarray(self.K_spin_free, dtype=complex)
        Hso = np.asarray(self.H_soc, dtype=complex)
        Kso = np.asarray(self.K_soc, dtype=complex)
        tolerance = float(tolerance)
        if not np.isfinite(tolerance) or tolerance <= 0.0:
            raise ValueError("component tolerance must be finite and positive.")
        if q.ndim != 1 or len(q) < 1:
            raise ValueError("SOC components require a nonempty coordinate vector.")
        if H0.ndim != 2 or H0.shape[0] < 1 or H0.shape[0] != H0.shape[1]:
            raise ValueError("spin-free SOC Hamiltonian component must be square.")
        nstate = H0.shape[0]
        if Hso.shape != (nstate, nstate):
            raise ValueError("spin-free and SOC Hamiltonian dimensions differ.")
        expected_derivative_shape = (len(q), nstate, nstate)
        if K0.shape != expected_derivative_shape or Kso.shape != expected_derivative_shape:
            raise ValueError("SOC derivatives have incompatible coordinate/state dimensions.")
        if not all(np.all(np.isfinite(item)) for item in (q, H0, K0, Hso, Kso)):
            raise ValueError("analytic SOC components contain non-finite data.")
        for name, matrix in (
            ("H_spin_free", H0),
            ("H_soc", Hso),
            ("H_total", H0 + Hso),
            *[(f"K_spin_free[{a}]", K0[a]) for a in range(len(q))],
            *[(f"K_soc[{a}]", Kso[a]) for a in range(len(q))],
            *[(f"K_total[{a}]", K0[a] + Kso[a]) for a in range(len(q))],
        ):
            if hermiticity_residual_v213(matrix) > tolerance:
                raise ValueError(f"{name} is not Hermitian.")
        return self


@dataclass(frozen=True)
class SingletTripletSOCConfigV220:
    mass: float = 900.0
    singlet_force_constant: float = 0.012
    triplet_force_constant: float = 0.010
    singlet_center: float = -0.55
    triplet_center: float = 0.55
    singlet_offset: float = 0.0
    triplet_offset: float = 0.004
    soc_scale: float = 0.003
    lambda_real_intercept: float = 0.75
    lambda_real_gradient: float = 0.10
    lambda_imag_intercept: float = -0.28
    lambda_imag_gradient: float = 0.06
    lambda_zero_intercept: float = 0.42
    lambda_zero_gradient: float = -0.04
    soc_enabled: bool = True

    def validate(self):
        _finite_real("mass", self.mass, positive=True)
        _finite_real("singlet_force_constant", self.singlet_force_constant, positive=True)
        _finite_real("triplet_force_constant", self.triplet_force_constant, positive=True)
        for name, value in asdict(self).items():
            if name not in {"mass", "singlet_force_constant", "triplet_force_constant", "soc_enabled"}:
                _finite_real(name, value)
        if not isinstance(self.soc_enabled, (bool, np.bool_)):
            raise ValueError("soc_enabled must be Boolean.")
        if not self.soc_enabled and float(self.soc_scale) != 0.0:
            raise ValueError("disabled SOC requires soc_scale=0 exactly.")
        return self

    def model_space(self, representation="fixed_spin_diabatic"):
        return ElectronicModelSpaceV213(
            name="v0.22.0 complete singlet-triplet analytic space",
            representation=representation,
            states=(
                ElectronicStateDescriptorV213("S", "S", 1, "M=0", 0),
                ElectronicStateDescriptorV213("T(M=-1)", "T", 3, "M=-1", 0),
                ElectronicStateDescriptorV213("T(M=0)", "T", 3, "M=0", 0),
                ElectronicStateDescriptorV213("T(M=+1)", "T", 3, "M=+1", 0),
            ),
            complete_multiplets=True,
        ).validate()

    def provenance(self, representation="fixed_spin_diabatic"):
        self.validate()
        symmetry = SOCSymmetryContractV221(
            "even",
            singlet_triplet_time_reversal_matrix_v220(),
            singlet_triplet_projectors_v220(),
        )
        return ElectronicOperatorProvenanceV213(
            model_name="v0.22 analytic singlet-triplet SOC",
            model_version="2",
            model_space=self.model_space(representation),
            spin_free_method="analytic displaced harmonic singlet/triplet surfaces",
            soc_enabled=bool(self.soc_enabled),
            soc_method=(
                "analytic time-reversal-invariant linear singlet-triplet SOC"
                if self.soc_enabled
                else "none"
            ),
            scalar_relativistic_method="none",
            derivative_method="analytic physical H_spin_free and H_SOC derivatives",
            parameters={
                **asdict(self),
                **symmetry.as_provenance_parameters(),
                "basis_order": ["S", "T(M=-1)", "T(M=0)", "T(M=+1)"],
                "time_reversal_convention": "Theta=J_ST K; J_ST J_ST*=+I",
                "soc_phase_convention": "(lambda, i*mu, lambda*)",
                "soc_signal_expected": bool(
                    self.soc_enabled and float(self.soc_scale) != 0.0
                ),
            },
        ).validate()


class AnalyticSingletTripletSOCProviderV220:
    def __init__(self, config=SingletTripletSOCConfigV220()):
        self.config = config.validate()
        self.provenance = self.config.provenance()
        self.calls = 0

    @property
    def time_reversal_matrix(self):
        return singlet_triplet_time_reversal_matrix_v220()

    @property
    def projectors(self):
        return singlet_triplet_projectors_v220()

    @property
    def soc_symmetry_contract(self):
        return SOCSymmetryContractV221(
            "even", self.time_reversal_matrix, self.projectors
        )

    def components(self, q):
        q = np.asarray(q, dtype=float)
        if q.shape != (1,) or not np.all(np.isfinite(q)):
            raise ValueError("singlet-triplet geometry must be one finite coordinate.")
        x = float(q[0])
        cfg = self.config
        Vs = 0.5 * cfg.singlet_force_constant * (x - cfg.singlet_center) ** 2 + cfg.singlet_offset
        Vt = 0.5 * cfg.triplet_force_constant * (x - cfg.triplet_center) ** 2 + cfg.triplet_offset
        dVs = cfg.singlet_force_constant * (x - cfg.singlet_center)
        dVt = cfg.triplet_force_constant * (x - cfg.triplet_center)
        H0 = np.diag([Vs, Vt, Vt, Vt]).astype(complex)
        K0 = np.diag([dVs, dVt, dVt, dVt]).astype(complex)[None, :, :]

        scale = float(cfg.soc_scale) if cfg.soc_enabled else 0.0
        lam = scale * (
            cfg.lambda_real_intercept
            + cfg.lambda_real_gradient * x
            + 1j * (cfg.lambda_imag_intercept + cfg.lambda_imag_gradient * x)
        )
        dlam = scale * (cfg.lambda_real_gradient + 1j * cfg.lambda_imag_gradient)
        mu = scale * (cfg.lambda_zero_intercept + cfg.lambda_zero_gradient * x)
        dmu = scale * cfg.lambda_zero_gradient
        Hso = np.zeros((4, 4), dtype=complex)
        Kso = np.zeros((1, 4, 4), dtype=complex)
        couplings = np.asarray([lam, 1j * mu, np.conj(lam)])
        derivatives = np.asarray([dlam, 1j * dmu, np.conj(dlam)])
        Hso[0, 1:] = couplings
        Hso[1:, 0] = np.conj(couplings)
        Kso[0, 0, 1:] = derivatives
        Kso[0, 1:, 0] = np.conj(derivatives)
        return SOCOperatorComponentsV220(q.copy(), H0, K0, Hso, Kso).validate()

    def evaluate_snapshot(self, q):
        self.calls += 1
        components = self.components(q)
        point = compose_electronic_operator_v213(
            q=components.q,
            H_spin_free=components.H_spin_free,
            dH_spin_free_dq=components.K_spin_free,
            H_soc=components.H_soc,
            dH_soc_dq=components.K_soc,
            connection_q=np.zeros((1, 4, 4), dtype=complex),
            mass_matrix_q_au=np.asarray([[self.config.mass]]),
            provenance=self.provenance,
        )
        return ElectronicOperatorSnapshotV21(
            point=point,
            state_vectors=np.eye(4, dtype=complex),
            metadata={
                "provider": "AnalyticSingletTripletSOCProviderV220",
                "physical_soc": bool(self.config.soc_enabled),
                "electron_parity": "even",
                "provenance_fingerprint": self.provenance.fingerprint(),
            },
        ).validate()

    def evaluate(self, q):
        return self.evaluate_snapshot(q).point

    def snapshot_overlap(self, left, right):
        return np.eye(4, dtype=complex)

    def diagnostics_dict(self):
        return {
            "provider": "AnalyticSingletTripletSOCProviderV220",
            "calls": int(self.calls),
            "physical_soc": bool(self.config.soc_enabled),
            "provenance_fingerprint": self.provenance.fingerprint(),
        }


@dataclass(frozen=True)
class DoubletSOCConfigV220:
    mass: float = 950.0
    root_1_force_constant: float = 0.011
    root_2_force_constant: float = 0.013
    root_1_center: float = -0.45
    root_2_center: float = 0.50
    root_1_offset: float = 0.0
    root_2_offset: float = 0.005
    soc_scale: float = 0.0025
    a_real_intercept: float = 0.68
    a_real_gradient: float = 0.07
    a_imag_intercept: float = -0.31
    a_imag_gradient: float = 0.05
    b_real_intercept: float = 0.38
    b_real_gradient: float = -0.06
    b_imag_intercept: float = 0.22
    b_imag_gradient: float = 0.04
    soc_enabled: bool = True

    def validate(self):
        _finite_real("mass", self.mass, positive=True)
        _finite_real("root_1_force_constant", self.root_1_force_constant, positive=True)
        _finite_real("root_2_force_constant", self.root_2_force_constant, positive=True)
        for name, value in asdict(self).items():
            if name not in {"mass", "root_1_force_constant", "root_2_force_constant", "soc_enabled"}:
                _finite_real(name, value)
        if not isinstance(self.soc_enabled, (bool, np.bool_)):
            raise ValueError("soc_enabled must be Boolean.")
        if not self.soc_enabled and float(self.soc_scale) != 0.0:
            raise ValueError("disabled SOC requires soc_scale=0 exactly.")
        return self

    def model_space(self, representation="fixed_spin_diabatic"):
        return ElectronicModelSpaceV213(
            name="v0.22.0 two-complete-doublet analytic space",
            representation=representation,
            states=(
                ElectronicStateDescriptorV213("D1(+1/2)", "D1", 2, "M=+1/2", 0),
                ElectronicStateDescriptorV213("D1(-1/2)", "D1", 2, "M=-1/2", 0),
                ElectronicStateDescriptorV213("D2(+1/2)", "D2", 2, "M=+1/2", 0),
                ElectronicStateDescriptorV213("D2(-1/2)", "D2", 2, "M=-1/2", 0),
            ),
            complete_multiplets=True,
        ).validate()

    def provenance(self, representation="fixed_spin_diabatic"):
        self.validate()
        symmetry = SOCSymmetryContractV221(
            "odd",
            kramers_time_reversal_matrix_v220(),
            doublet_root_projectors_v220(),
        )
        return ElectronicOperatorProvenanceV213(
            model_name="v0.22 analytic two-Kramers-doublet SOC",
            model_version="2",
            model_space=self.model_space(representation),
            spin_free_method="analytic displaced harmonic doublet-root surfaces",
            soc_enabled=bool(self.soc_enabled),
            soc_method=(
                "analytic time-reversal-invariant quaternionic doublet SOC"
                if self.soc_enabled
                else "none"
            ),
            scalar_relativistic_method="none",
            derivative_method="analytic physical H_spin_free and H_SOC derivatives",
            parameters={
                **asdict(self),
                **symmetry.as_provenance_parameters(),
                "basis_order": ["D1(+1/2)", "D1(-1/2)", "D2(+1/2)", "D2(-1/2)"],
                "time_reversal_convention": "Theta=J_D K; J_D J_D*=-I",
                "soc_block_convention": "B=[[a,b],[-b*,a*]]",
                "soc_signal_expected": bool(
                    self.soc_enabled and float(self.soc_scale) != 0.0
                ),
            },
        ).validate()


def _quaternion_block_v220(a, b):
    return np.asarray([[a, b], [-np.conj(b), np.conj(a)]], dtype=complex)


class AnalyticDoubletSOCProviderV220:
    def __init__(self, config=DoubletSOCConfigV220()):
        self.config = config.validate()
        self.provenance = self.config.provenance()
        self.calls = 0

    @property
    def time_reversal_matrix(self):
        return kramers_time_reversal_matrix_v220()

    @property
    def projectors(self):
        return doublet_root_projectors_v220()

    @property
    def soc_symmetry_contract(self):
        return SOCSymmetryContractV221(
            "odd", self.time_reversal_matrix, self.projectors
        )

    def components(self, q):
        q = np.asarray(q, dtype=float)
        if q.shape != (1,) or not np.all(np.isfinite(q)):
            raise ValueError("doublet geometry must be one finite coordinate.")
        x = float(q[0])
        cfg = self.config
        E1 = 0.5 * cfg.root_1_force_constant * (x - cfg.root_1_center) ** 2 + cfg.root_1_offset
        E2 = 0.5 * cfg.root_2_force_constant * (x - cfg.root_2_center) ** 2 + cfg.root_2_offset
        dE1 = cfg.root_1_force_constant * (x - cfg.root_1_center)
        dE2 = cfg.root_2_force_constant * (x - cfg.root_2_center)
        H0 = np.diag([E1, E1, E2, E2]).astype(complex)
        K0 = np.diag([dE1, dE1, dE2, dE2]).astype(complex)[None, :, :]

        scale = float(cfg.soc_scale) if cfg.soc_enabled else 0.0
        a = scale * (
            cfg.a_real_intercept + cfg.a_real_gradient * x
            + 1j * (cfg.a_imag_intercept + cfg.a_imag_gradient * x)
        )
        b = scale * (
            cfg.b_real_intercept + cfg.b_real_gradient * x
            + 1j * (cfg.b_imag_intercept + cfg.b_imag_gradient * x)
        )
        da = scale * (cfg.a_real_gradient + 1j * cfg.a_imag_gradient)
        db = scale * (cfg.b_real_gradient + 1j * cfg.b_imag_gradient)
        block = _quaternion_block_v220(a, b)
        dblock = _quaternion_block_v220(da, db)
        Hso = np.zeros((4, 4), dtype=complex)
        Kso = np.zeros((1, 4, 4), dtype=complex)
        Hso[:2, 2:] = block
        Hso[2:, :2] = block.conj().T
        Kso[0, :2, 2:] = dblock
        Kso[0, 2:, :2] = dblock.conj().T
        return SOCOperatorComponentsV220(q.copy(), H0, K0, Hso, Kso).validate()

    def evaluate_snapshot(self, q):
        self.calls += 1
        components = self.components(q)
        point = compose_electronic_operator_v213(
            q=components.q,
            H_spin_free=components.H_spin_free,
            dH_spin_free_dq=components.K_spin_free,
            H_soc=components.H_soc,
            dH_soc_dq=components.K_soc,
            connection_q=np.zeros((1, 4, 4), dtype=complex),
            mass_matrix_q_au=np.asarray([[self.config.mass]]),
            provenance=self.provenance,
        )
        return ElectronicOperatorSnapshotV21(
            point=point,
            state_vectors=np.eye(4, dtype=complex),
            metadata={
                "provider": "AnalyticDoubletSOCProviderV220",
                "physical_soc": bool(self.config.soc_enabled),
                "electron_parity": "odd",
                "provenance_fingerprint": self.provenance.fingerprint(),
            },
        ).validate()

    def evaluate(self, q):
        return self.evaluate_snapshot(q).point

    def snapshot_overlap(self, left, right):
        return np.eye(4, dtype=complex)

    def diagnostics_dict(self):
        return {
            "provider": "AnalyticDoubletSOCProviderV220",
            "calls": int(self.calls),
            "physical_soc": bool(self.config.soc_enabled),
            "provenance_fingerprint": self.provenance.fingerprint(),
        }
