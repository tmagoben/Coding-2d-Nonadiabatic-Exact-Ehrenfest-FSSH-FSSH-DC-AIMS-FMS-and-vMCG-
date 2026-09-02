import numpy as np

from .electronic_structure import (
    ElectronicStructurePoint,
    project_cartesian_vector_to_coordinate,
)
from .pyscf_nac_convention_v232 import (
    PYSCF_NAC_EMPIRICAL_MAPPING_V232,
    PYSCF_NAC_UPSTREAM_DOCUMENTATION_V232,
    pyscf_state_tuple_for_internal_dij_v232,
    require_exact_pyscf_version_v232,
)
from .nac_compatibility_v233 import corrected_pyscf_nac_convention_v233


class PySCFStateAveragedCASSCFProvider:
    """Optional 1D provider based on PySCF SA-CASSCF gradients and NACs.

    geometry_builder(q) must return
        atoms, tangent

    atoms:
        PySCF atom specification, preferably a list such as
        [("H",(0,0,-R/2)), ("H",(0,0,R/2))]

    tangent:
        ndarray (natm,3) giving dR_cartesian/dq in the same coordinate units.

    Conventions
    -----------
    This provider stores
        nac_q[i,j] = <phi_i | d/dq phi_j>.

    PySCF's SA-CASSCF NAC API documents ``state=(ket,bra)`` as returning
    ``<bra|d(ket)/dR>``.  PySCF 2.13.1 central differences of phase-aligned
    many-electron overlaps empirically require ``state=(i,j)`` for the internal
    ``d[i,j]`` convention.  v0.23.2 uses that certified production mapping.
    """
    def __init__(
        self,
        geometry_builder,
        basis,
        ncas,
        nelecas,
        nstates=2,
        charge=0,
        spin=0,
        scf_conv_tol=1e-10,
        mc_conv_tol=1e-9,
        mc_conv_tol_grad=1e-5,
        use_etfs=False,
        verbose=0,
    ):
        self.geometry_builder=geometry_builder
        self.basis=basis
        self.ncas=int(ncas)
        self.nelecas=nelecas
        self.nstates=int(nstates)
        self.charge=int(charge)
        self.spin=int(spin)
        self.scf_conv_tol=float(scf_conv_tol)
        self.mc_conv_tol=float(mc_conv_tol)
        self.mc_conv_tol_grad=float(mc_conv_tol_grad)
        self.use_etfs=bool(use_etfs)
        self.verbose=int(verbose)

    @staticmethod
    def _imports():
        try:
            import pyscf
            from pyscf import gto, scf, mcscf
        except ImportError as exc:
            raise ImportError(
                "PySCF is optional. Install this project with `pip install -e '.[pyscf]'`."
            ) from exc
        require_exact_pyscf_version_v232(pyscf)
        return gto,scf,mcscf

    def evaluate(self,q):
        gto,scf,mcscf=self._imports()
        nac_identity=corrected_pyscf_nac_convention_v233(use_etfs=self.use_etfs)
        atoms,tangent=self.geometry_builder(float(q))
        tangent=np.asarray(tangent,float)

        mol=gto.M(
            atom=atoms,
            basis=self.basis,
            charge=self.charge,
            spin=self.spin,
            unit="Bohr",
            verbose=self.verbose,
        )

        mf=scf.RHF(mol)
        mf.conv_tol=self.scf_conv_tol
        mf.kernel()
        if not mf.converged:
            raise RuntimeError(f"RHF did not converge at q={q}.")

        weights=np.ones(self.nstates)/self.nstates
        mc=mcscf.CASSCF(mf,self.ncas,self.nelecas).state_average(weights)
        mc.conv_tol=self.mc_conv_tol
        mc.conv_tol_grad=self.mc_conv_tol_grad
        mc.kernel()
        if not getattr(mc,"converged",False):
            raise RuntimeError(f"SA-CASSCF did not converge at q={q}.")

        energies=np.asarray(mc.e_states,float)
        gradients=np.zeros(self.nstates,float)
        for state in range(self.nstates):
            grad_cart=np.asarray(mc.nuc_grad_method(state=state).kernel(),float)
            gradients[state]=project_cartesian_vector_to_coordinate(grad_cart,tangent)

        nac=np.zeros((self.nstates,self.nstates),float)
        nac_method=mc.nac_method()
        for i in range(self.nstates):
            for j in range(i+1,self.nstates):
                cart=np.asarray(
                    nac_method.kernel(
                        state=pyscf_state_tuple_for_internal_dij_v232(i,j),
                        use_etfs=self.use_etfs,
                    ),
                    float,
                )
                value=project_cartesian_vector_to_coordinate(cart,tangent)
                nac[i,j]=value
                nac[j,i]=-value

        return ElectronicStructurePoint(
            q=float(q),
            energies=energies,
            gradients_q=gradients,
            nac_q=nac,
            metadata={
                "provider":"pyscf_sa_casscf",
                "basis":self.basis,
                "ncas":self.ncas,
                "nelecas":self.nelecas,
                "nstates":self.nstates,
                "weights":weights.tolist(),
                "scf_converged":bool(mf.converged),
                "mc_converged":bool(mc.converged),
                "use_etfs":self.use_etfs,
                "nac_convention":"nac_q[i,j]=<phi_i|d/dq phi_j>",
                "pyscf_state_tuple":"state=(ket,bra)",
                "pyscf_request_for_internal_dij":"state=(i,j)",
                "pyscf_nac_upstream_documentation":(
                    PYSCF_NAC_UPSTREAM_DOCUMENTATION_V232
                ),
                "pyscf_nac_empirical_mapping":PYSCF_NAC_EMPIRICAL_MAPPING_V232,
                "pyscf_nac_mapping_certification":(
                    "v0.23.2 phase-aligned many-electron overlap central difference"
                ),
                "v233_nac_convention":nac_identity.as_dict(),
                "v233_nac_convention_fingerprint":nac_identity.fingerprint(),
            },
        ).validate()
