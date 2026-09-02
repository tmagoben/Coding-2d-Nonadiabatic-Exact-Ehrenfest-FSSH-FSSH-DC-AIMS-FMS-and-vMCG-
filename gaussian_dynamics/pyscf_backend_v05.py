from dataclasses import dataclass
import json
import hashlib
import numpy as np

from .molecular_backend import (
    MolecularGeometry,
    CartesianElectronicStructurePoint,
)
from .pyscf_nac_convention_v232 import (
    PYSCF_NAC_EMPIRICAL_MAPPING_V232,
    PYSCF_NAC_UPSTREAM_DOCUMENTATION_V232,
    pyscf_state_tuple_for_internal_dij_v232,
    require_exact_pyscf_version_v232,
)
from .nac_compatibility_v233 import corrected_pyscf_nac_convention_v233


@dataclass(frozen=True)
class PySCFSACASSCFConfig:
    basis: object
    ncas: int
    nelecas: object
    nstates: int = 2
    weights: tuple | None = None

    charge: int = 0
    spin: int = 0
    symmetry: bool = False
    scf_reference: str = "RHF"

    scf_conv_tol: float = 1e-10
    scf_max_cycle: int = 100

    mc_conv_tol: float = 1e-9
    mc_conv_tol_grad: float = 1e-5
    mc_max_cycle_macro: int = 50

    use_etfs: bool = False
    compute_scaled_nac: bool = False
    warm_start_mo: bool = False
    isotope_avg_masses: bool = True

    verbose: int = 0
    max_memory_mb: int = 2000

    def normalized_weights(self):
        if self.nstates < 2:
            raise ValueError("SA-CASSCF backend requires at least two states.")
        if self.ncas <= 0:
            raise ValueError("ncas must be positive.")

        if self.weights is None:
            w = np.ones(self.nstates, dtype=float) / self.nstates
        else:
            w = np.asarray(self.weights, dtype=float)

        if w.shape != (self.nstates,):
            raise ValueError("weights must have length nstates.")
        if np.any(w < 0.0) or not np.isclose(np.sum(w), 1.0):
            raise ValueError("weights must be nonnegative and sum to one.")
        return w

    def fingerprint(self):
        payload = {
            "basis": repr(self.basis),
            "ncas": self.ncas,
            "nelecas": repr(self.nelecas),
            "nstates": self.nstates,
            "weights": self.normalized_weights().tolist(),
            "charge": self.charge,
            "spin": self.spin,
            "symmetry": self.symmetry,
            "scf_reference": self.scf_reference,
            "scf_conv_tol": self.scf_conv_tol,
            "scf_max_cycle": self.scf_max_cycle,
            "mc_conv_tol": self.mc_conv_tol,
            "mc_conv_tol_grad": self.mc_conv_tol_grad,
            "mc_max_cycle_macro": self.mc_max_cycle_macro,
            "use_etfs": self.use_etfs,
            "compute_scaled_nac": self.compute_scaled_nac,
            "isotope_avg_masses": self.isotope_avg_masses,
        }
        raw = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
        return hashlib.sha256(raw).hexdigest()


class PySCFSACASSCFBackend:
    """Explicit PySCF SA-CASSCF energies, gradients, and NAC backend.

    Internal contract
    -----------------
    nac_cart[i,j] = <Phi_i | grad_R Phi_j>.

    PySCF 2.13.1 empirical contract
    -------------------------------
    Although the upstream docstring describes ``state=(ket,bra)`` as returning
    ``<bra|grad_R ket>``, central differences of phase-aligned many-electron
    overlaps show that internal ``nac_cart[i,j]=<i|grad_R j>`` is obtained with
    ``state=(i,j)``.  v0.23.2 freezes that runtime-certified mapping.
    """

    def __init__(self, config: PySCFSACASSCFConfig):
        self.config = config
        self.config.normalized_weights()
        if self.config.scf_reference.upper() not in {"RHF", "ROHF"}:
            raise ValueError("scf_reference must be 'RHF' or 'ROHF'.")
        self._last_mo = None

    @staticmethod
    def _imports():
        try:
            import pyscf
            from pyscf import gto, scf, mcscf
        except ImportError as exc:
            raise ImportError(
                "PySCF is required for PySCFSACASSCFBackend. "
                "Install with `pip install -e '.[pyscf]'`."
            ) from exc
        require_exact_pyscf_version_v232(pyscf)
        return pyscf, gto, scf, mcscf

    def _build_molecule(self, geometry, gto):
        atoms = [
            (symbol, tuple(float(v) for v in coord))
            for symbol, coord in zip(geometry.symbols, geometry.coords_bohr)
        ]

        return gto.M(
            atom=atoms,
            basis=self.config.basis,
            charge=self.config.charge,
            spin=self.config.spin,
            unit="Bohr",
            symmetry=self.config.symmetry,
            verbose=self.config.verbose,
            max_memory=self.config.max_memory_mb,
        )

    def evaluate(self, geometry):
        geometry = MolecularGeometry(geometry.symbols, geometry.coords_bohr)
        pyscf, gto, scf, mcscf = self._imports()

        mol = self._build_molecule(geometry, gto)

        ref = self.config.scf_reference.upper()
        mf_cls = scf.RHF if ref == "RHF" else scf.ROHF
        mf = mf_cls(mol)
        mf.conv_tol = self.config.scf_conv_tol
        mf.max_cycle = self.config.scf_max_cycle
        mf.kernel()

        if not bool(getattr(mf, "converged", False)):
            raise RuntimeError("PySCF SCF did not converge.")

        weights = self.config.normalized_weights()

        mc = mcscf.CASSCF(mf, self.config.ncas, self.config.nelecas)
        mc = mc.state_average_(weights)
        mc.conv_tol = self.config.mc_conv_tol
        mc.conv_tol_grad = self.config.mc_conv_tol_grad
        mc.max_cycle_macro = self.config.mc_max_cycle_macro

        if self.config.warm_start_mo and self._last_mo is not None:
            mc.kernel(self._last_mo)
        else:
            mc.kernel()

        if not bool(getattr(mc, "converged", False)):
            raise RuntimeError("PySCF SA-CASSCF did not converge.")

        self._last_mo = np.asarray(mc.mo_coeff).copy()

        energies = np.asarray(mc.e_states, dtype=float)
        if energies.shape != (self.config.nstates,):
            raise RuntimeError(
                f"Expected {self.config.nstates} SA-CASSCF state energies, "
                f"received shape {energies.shape}."
            )

        gradients = np.zeros(
            (self.config.nstates, geometry.natom, 3),
            dtype=float,
        )

        for state in range(self.config.nstates):
            grad_method = mc.nuc_grad_method(state=state)
            gradients[state] = np.asarray(grad_method.kernel(), dtype=float)

        nac = np.zeros(
            (self.config.nstates, self.config.nstates, geometry.natom, 3),
            dtype=float,
        )

        scaled = None
        if self.config.compute_scaled_nac:
            scaled = np.zeros_like(nac)

        nac_method = mc.nac_method()

        for i in range(self.config.nstates):
            for j in range(i + 1, self.config.nstates):
                # PySCF 2.13.1 empirical mapping, certified against central
                # differences of phase-aligned many-electron overlaps.
                dij = np.asarray(
                    nac_method.kernel(
                        state=pyscf_state_tuple_for_internal_dij_v232(i, j),
                        use_etfs=self.config.use_etfs,
                        mult_ediff=False,
                    ),
                    dtype=float,
                )

                nac[i, j] = dij
                nac[j, i] = -dij

                if scaled is not None:
                    sij = np.asarray(
                        nac_method.kernel(
                            state=pyscf_state_tuple_for_internal_dij_v232(i, j),
                            use_etfs=self.config.use_etfs,
                            mult_ediff=True,
                        ),
                        dtype=float,
                    )
                    # PySCF's mult_ediff quantity is symmetric under swapping the
                    # two state arguments in its own convention.
                    scaled[i, j] = sij
                    scaled[j, i] = sij

        masses = np.asarray(
            mol.atom_mass_list(isotope_avg=self.config.isotope_avg_masses),
            dtype=float,
        )

        nac_identity = corrected_pyscf_nac_convention_v233(
            use_etfs=self.config.use_etfs
        )
        metadata = {
            "backend": "PySCF SA-CASSCF",
            "pyscf_version": getattr(pyscf, "__version__", "unknown"),
            "basis": repr(self.config.basis),
            "charge": self.config.charge,
            "spin": self.config.spin,
            "symmetry": self.config.symmetry,
            "scf_reference": ref,
            "ncas": self.config.ncas,
            "nelecas": repr(self.config.nelecas),
            "nstates": self.config.nstates,
            "state_average_weights": weights.tolist(),
            "scf_conv_tol": self.config.scf_conv_tol,
            "scf_max_cycle": self.config.scf_max_cycle,
            "mc_conv_tol": self.config.mc_conv_tol,
            "mc_conv_tol_grad": self.config.mc_conv_tol_grad,
            "mc_max_cycle_macro": self.config.mc_max_cycle_macro,
            "scf_converged": True,
            "mc_converged": True,
            "use_etfs": self.config.use_etfs,
            "compute_scaled_nac": self.config.compute_scaled_nac,
            "nac_internal_convention": "d[i,j]=<Phi_i|grad_R Phi_j>",
            "pyscf_nac_state_tuple": "(ket,bra)",
            "pyscf_request_for_internal_dij": "state=(i,j)",
            "pyscf_nac_upstream_documentation": (
                PYSCF_NAC_UPSTREAM_DOCUMENTATION_V232
            ),
            "pyscf_nac_empirical_mapping": PYSCF_NAC_EMPIRICAL_MAPPING_V232,
            "pyscf_nac_mapping_certification": (
                "v0.23.2 phase-aligned many-electron overlap central difference"
            ),
            "v233_nac_convention": nac_identity.as_dict(),
            "v233_nac_convention_fingerprint": nac_identity.fingerprint(),
            "dynamics_mult_ediff": False,
            "warm_start_mo": self.config.warm_start_mo,
            "state_tracking": "energy/root ordering only; no many-electron overlap tracking",
            "config_fingerprint": self.config.fingerprint(),
        }

        if scaled is not None:
            metadata["scaled_nac_semantics"] = (
                "Raw PySCF mult_ediff=True value requested with state=(i,j); "
                "stored for diagnostics only and not used as the dynamics derivative coupling."
            )

        return CartesianElectronicStructurePoint(
            geometry=geometry,
            energies=energies,
            gradients_cart=gradients,
            nac_cart=nac,
            masses_amu=masses,
            scaled_nac_cart=scaled,
            metadata=metadata,
        ).validate()
