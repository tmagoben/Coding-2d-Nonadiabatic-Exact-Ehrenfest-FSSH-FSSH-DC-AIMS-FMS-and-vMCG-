import warnings
import numpy as np

from .molecular_backend import (
    MolecularGeometry,
    CartesianElectronicStructurePoint,
)
from .pyscf_backend_v05 import PySCFSACASSCFConfig
from .pyscf_wavefunction_overlap import (
    CASSCFWavefunctionSnapshot,
    casscf_state_overlap_matrix,
)
from .pyscf_nac_convention_v232 import (
    PYSCF_NAC_EMPIRICAL_MAPPING_V232,
    PYSCF_NAC_UPSTREAM_DOCUMENTATION_V232,
    pyscf_state_tuple_for_internal_dij_v232,
    require_exact_pyscf_version_v232,
)
from .nac_compatibility_v233 import corrected_pyscf_nac_convention_v233
from .state_tracking import (
    maximum_overlap_assignment,
    transform_state_properties,
    energy_degeneracy_clusters,
)


class PySCFTrackedSACASSCFBackend:
    """State-overlap-tracked PySCF SA-CASSCF backend for a sequential geometry path.

    This object is intentionally *stateful*.  Each call is tracked against the
    immediately preceding accepted geometry.

    It is appropriate for:
      - geometry scans,
      - one TBF trajectory history,
      - sequential direct-dynamics center propagation.

    It should not be shared blindly among unrelated/asynchronously sampled Gaussian
    pair centroids because "previous geometry" would then depend on call order.
    """

    def __init__(
        self,
        config: PySCFSACASSCFConfig,
        minimum_overlap=0.50,
        minimum_score_margin=0.05,
        energy_degeneracy_tolerance=1e-4,
        ambiguity_policy="raise",
        overlap_engine=None,
    ):
        self.config = config
        self.config.normalized_weights()

        if self.config.scf_reference.upper() not in {"RHF", "ROHF"}:
            raise ValueError("scf_reference must be 'RHF' or 'ROHF'.")

        if ambiguity_policy not in {"raise", "warn", "accept"}:
            raise ValueError("ambiguity_policy must be 'raise', 'warn', or 'accept'.")

        self.minimum_overlap = float(minimum_overlap)
        self.minimum_score_margin = float(minimum_score_margin)
        self.energy_degeneracy_tolerance = float(energy_degeneracy_tolerance)
        self.ambiguity_policy = ambiguity_policy
        self.overlap_engine = overlap_engine or casscf_state_overlap_matrix

        self._previous_snapshot = None
        self._last_mo = None
        self.step_index = 0
        self.history = []

    @staticmethod
    def _imports():
        try:
            import pyscf
            from pyscf import gto, scf, mcscf
        except ImportError as exc:
            raise ImportError(
                "PySCF is required for PySCFTrackedSACASSCFBackend. "
                "Install with `pip install -e '.[pyscf]'`."
            ) from exc
        require_exact_pyscf_version_v232(pyscf)
        return pyscf, gto, scf, mcscf

    def reset_tracking(self, reset_orbitals=True):
        self._previous_snapshot = None
        self.step_index = 0
        self.history = []
        if reset_orbitals:
            self._last_mo = None

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

    def _run_raw(self, geometry):
        """Run one raw energy-ordered PySCF SA-CASSCF calculation + snapshot."""
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
                f"Expected {self.config.nstates} state energies, got {energies.shape}."
            )

        ci_roots = tuple(np.asarray(c).copy() for c in mc.ci)
        if len(ci_roots) != self.config.nstates:
            raise RuntimeError(
                f"Expected {self.config.nstates} CI roots, got {len(ci_roots)}."
            )

        if isinstance(mc.ncore, tuple):
            raise NotImplementedError(
                "v0.6 many-electron overlap tracking currently supports restricted "
                "RHF/ROHF CASSCF with a scalar ncore, not UHF-CASSCF."
            )

        gradients = np.zeros(
            (self.config.nstates, geometry.natom, 3),
            dtype=float,
        )
        for state in range(self.config.nstates):
            gradients[state] = np.asarray(
                mc.nuc_grad_method(state=state).kernel(),
                dtype=float,
            )

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
            "backend": "PySCF tracked SA-CASSCF v0.6",
            "pyscf_version": getattr(pyscf, "__version__", "unknown"),
            "basis": repr(self.config.basis),
            "charge": self.config.charge,
            "spin": self.config.spin,
            "scf_reference": ref,
            "ncas": int(mc.ncas),
            "ncore": int(mc.ncore),
            "nelecas": [int(x) for x in mc.nelecas],
            "nstates": self.config.nstates,
            "state_average_weights": weights.tolist(),
            "use_etfs": self.config.use_etfs,
            "dynamics_mult_ediff": False,
            "nac_internal_convention": "d[i,j]=<Phi_i|grad_R Phi_j>",
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
            "warm_start_mo": self.config.warm_start_mo,
            "raw_state_order": list(range(self.config.nstates)),
        }

        point = CartesianElectronicStructurePoint(
            geometry=geometry,
            energies=energies,
            gradients_cart=gradients,
            nac_cart=nac,
            masses_amu=masses,
            scaled_nac_cart=scaled,
            metadata=metadata,
        ).validate()

        snapshot = CASSCFWavefunctionSnapshot(
            mol=mol,
            mo_coeff=np.asarray(mc.mo_coeff).copy(),
            ci_roots=ci_roots,
            ncore=int(mc.ncore),
            ncas=int(mc.ncas),
            nelecas=tuple(int(x) for x in mc.nelecas),
            metadata={
                "pyscf_version": getattr(pyscf, "__version__", "unknown"),
                "step_index_raw": self.step_index,
            },
        )

        return point, snapshot

    def _handle_ambiguity(self, result):
        if not result.ambiguous:
            return

        message = (
            "Electronic state tracking is ambiguous: "
            + "; ".join(result.reasons)
        )

        if self.ambiguity_policy == "raise":
            raise RuntimeError(message)
        if self.ambiguity_policy == "warn":
            warnings.warn(message, RuntimeWarning, stacklevel=2)

    def evaluate_raw_with_snapshot(self, geometry):
        """Public nontracking PySCF point + CASSCF snapshot for graph dynamics."""
        return self._run_raw(geometry)

    def evaluate(self, geometry):
        raw_point, raw_snapshot = self._run_raw(geometry)

        raw_energies = raw_point.energies.copy()
        raw_clusters = energy_degeneracy_clusters(
            raw_energies,
            tolerance=self.energy_degeneracy_tolerance,
        )

        if self._previous_snapshot is None:
            # The first geometry defines tracked labels and phases.
            tracked_point = raw_point
            tracked_snapshot = raw_snapshot

            tracking_metadata = {
                "tracking_step": int(self.step_index),
                "tracking_reference": "initial geometry",
                "permutation_tracked_to_raw": list(range(self.config.nstates)),
                "phase_factors": [[1.0, 0.0]] * self.config.nstates,
                "assigned_overlap_magnitudes": [1.0] * self.config.nstates,
                "ambiguous": False,
                "raw_energy_degeneracy_clusters": [list(c) for c in raw_clusters],
            }
            overlap = np.eye(self.config.nstates, dtype=complex)

        else:
            overlap = np.asarray(
                self.overlap_engine(self._previous_snapshot, raw_snapshot),
                dtype=complex,
            )

            result = maximum_overlap_assignment(
                overlap,
                minimum_overlap=self.minimum_overlap,
                minimum_score_margin=self.minimum_score_margin,
                real_gauge=True,
            )

            self._handle_ambiguity(result)

            E, G, D = transform_state_properties(
                raw_point.energies,
                raw_point.gradients_cart,
                raw_point.nac_cart,
                result,
            )

            scaled = None
            if raw_point.scaled_nac_cart is not None:
                # Diagnostic scaled NAC: same root permutation and real-state signs.
                perm = result.permutation
                phase = np.real(result.phase_factors)
                scaled = raw_point.scaled_nac_cart[np.ix_(perm, perm)].copy()
                sign_matrix = phase[:, None] * phase[None, :]
                scaled *= sign_matrix[..., None, None]

            tracked_metadata = dict(raw_point.metadata)
            tracked_metadata.update(result.as_metadata())
            tracked_overlap = (
                overlap[:, result.permutation]
                * result.phase_factors[None, :]
            )

            singular_values = np.linalg.svd(tracked_overlap, compute_uv=False)
            unitarity_defect = np.linalg.norm(
                tracked_overlap.conj().T @ tracked_overlap
                - np.eye(tracked_overlap.shape[1]),
                ord="fro",
            )

            tracked_metadata.update({
                "tracking_step": int(self.step_index),
                "tracking_reference": "previous accepted geometry",
                "raw_energies": raw_point.energies.tolist(),
                "tracked_energies": np.asarray(E).tolist(),
                "raw_energy_degeneracy_clusters": [list(c) for c in raw_clusters],
                "state_overlap_singular_values": singular_values.tolist(),
                "state_overlap_unitarity_defect": float(unitarity_defect),
                "state_overlap_matrix_raw": [
                    [
                        [float(np.real(z)), float(np.imag(z))]
                        for z in row
                    ]
                    for row in overlap
                ],
                "state_overlap_matrix_tracked_gauge": [
                    [
                        [float(np.real(z)), float(np.imag(z))]
                        for z in row
                    ]
                    for row in tracked_overlap
                ],
            })

            tracked_point = CartesianElectronicStructurePoint(
                geometry=raw_point.geometry,
                energies=np.asarray(E, dtype=float),
                gradients_cart=np.asarray(G, dtype=float),
                nac_cart=np.asarray(D, dtype=float),
                masses_amu=raw_point.masses_amu.copy(),
                scaled_nac_cart=scaled,
                metadata=tracked_metadata,
            ).validate()

            tracked_snapshot = raw_snapshot.with_transformed_roots(
                result.permutation,
                result.phase_factors,
            )
            tracking_metadata = tracked_metadata

        self._previous_snapshot = tracked_snapshot

        record = {
            "step": int(self.step_index),
            "energies": tracked_point.energies.tolist(),
            "tracking": tracking_metadata,
        }
        self.history.append(record)
        self.step_index += 1

        return tracked_point

    @property
    def previous_snapshot(self):
        return self._previous_snapshot
