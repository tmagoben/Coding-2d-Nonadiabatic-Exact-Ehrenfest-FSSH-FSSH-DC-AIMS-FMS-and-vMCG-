from dataclasses import dataclass
import numpy as np

from .pyscf_nac_convention_v232 import require_exact_pyscf_version_v232


@dataclass
class CASSCFWavefunctionSnapshot:
    """Minimal restricted-CASSCF data required for cross-geometry state overlaps."""
    mol: object
    mo_coeff: np.ndarray
    ci_roots: tuple
    ncore: int
    ncas: int
    nelecas: tuple
    metadata: dict

    def __post_init__(self):
        self.mo_coeff = np.asarray(self.mo_coeff)
        self.ci_roots = tuple(np.asarray(c) for c in self.ci_roots)
        self.ncore = int(self.ncore)
        self.ncas = int(self.ncas)
        self.nelecas = tuple(int(x) for x in self.nelecas)

        if len(self.nelecas) != 2:
            raise ValueError("nelecas must be (nalpha_active, nbeta_active).")
        if self.ncore < 0 or self.ncas <= 0:
            raise ValueError("Invalid ncore/ncas.")
        if self.mo_coeff.ndim != 2:
            raise ValueError("Restricted CASSCF snapshot expects one MO coefficient matrix.")
        if self.mo_coeff.shape[1] < self.ncore + self.ncas:
            raise ValueError("MO coefficient matrix does not contain core + active orbitals.")

    @property
    def nroots(self):
        return len(self.ci_roots)

    @property
    def norb_correlated(self):
        return self.ncore + self.ncas

    @property
    def nelec_correlated(self):
        na, nb = self.nelecas
        return (self.ncore + na, self.ncore + nb)

    @property
    def correlated_mo_coeff(self):
        return self.mo_coeff[:, : self.ncore + self.ncas]

    def with_transformed_roots(self, permutation, phase_factors):
        permutation = np.asarray(permutation, dtype=int)
        phase_factors = np.asarray(phase_factors, dtype=complex)

        roots = []
        for i, raw_index in enumerate(permutation):
            roots.append(
                np.asarray(self.ci_roots[int(raw_index)], dtype=complex)
                * phase_factors[i]
            )

        meta = dict(self.metadata)
        meta["state_transform_permutation"] = permutation.tolist()
        meta["state_transform_phase"] = [
            [float(np.real(z)), float(np.imag(z))]
            for z in phase_factors
        ]

        return CASSCFWavefunctionSnapshot(
            mol=self.mol,
            mo_coeff=self.mo_coeff.copy(),
            ci_roots=tuple(roots),
            ncore=self.ncore,
            ncas=self.ncas,
            nelecas=self.nelecas,
            metadata=meta,
        )


def _imports():
    try:
        import pyscf
        from pyscf import gto, fci
        from pyscf.fci import cistring
    except ImportError as exc:
        raise ImportError(
            "PySCF is required for many-electron CASSCF state-overlap tracking. "
            "Install with `pip install -e '.[pyscf]'`."
        ) from exc

    require_exact_pyscf_version_v232(pyscf)

    return gto, fci, cistring


def embed_active_ci_with_doubly_occupied_core(
    ci_active,
    ncore,
    ncas,
    nelecas,
    cistring=None,
):
    """Embed a CAS CI vector into the full core+active FCI determinant space.

    The orbital ordering is

        [doubly occupied core | active].

    Core orbitals are occupied in every alpha and beta determinant. Virtual orbitals
    are omitted because their occupation is identically zero in the CASSCF wavefunction.

    This embedding lets `pyscf.fci.addons.overlap` include core-active cross overlaps
    between two geometries rather than factorizing the overlap into an approximate
    core determinant times an active-space CI overlap.
    """
    if cistring is None:
        _, _, cistring = _imports()

    ncore = int(ncore)
    ncas = int(ncas)
    neleca, nelecb = tuple(int(x) for x in nelecas)

    norb = ncore + ncas
    na_total = ncore + neleca
    nb_total = ncore + nelecb

    if norb >= 64:
        raise NotImplementedError(
            "The explicit determinant-string embedding currently requires "
            "ncore+ncas < 64."
        )

    ci_active = np.asarray(ci_active)

    active_a_strings = cistring.make_strings(range(ncas), neleca)
    active_b_strings = cistring.make_strings(range(ncas), nelecb)

    if ci_active.shape != (len(active_a_strings), len(active_b_strings)):
        raise ValueError(
            "Active CI vector shape is inconsistent with ncas/nelecas: "
            f"received {ci_active.shape}, expected "
            f"({len(active_a_strings)}, {len(active_b_strings)})."
        )

    core_mask = (1 << ncore) - 1 if ncore else 0

    full_a_strings = core_mask | (active_a_strings.astype(np.int64) << ncore)
    full_b_strings = core_mask | (active_b_strings.astype(np.int64) << ncore)

    addr_a = cistring.strs2addr(norb, na_total, full_a_strings)
    addr_b = cistring.strs2addr(norb, nb_total, full_b_strings)

    nstr_a = len(cistring.make_strings(range(norb), na_total))
    nstr_b = len(cistring.make_strings(range(norb), nb_total))

    full = np.zeros((nstr_a, nstr_b), dtype=ci_active.dtype)
    full[np.ix_(addr_a, addr_b)] = ci_active
    return full


def casscf_state_overlap_matrix(previous, current):
    """Compute cross-geometry many-electron overlaps for restricted CASSCF roots.

    O_ij = <Psi_i(previous) | Psi_j(current)>

    The calculation uses:
      1. PySCF cross-AO overlap integrals,
      2. the previous/current core+active MO overlap matrix,
      3. exact embedding of each CAS CI vector into the core+active FCI space,
      4. `pyscf.fci.addons.overlap(..., s=S_mo)`.

    Within the restricted CASSCF wavefunction represented by the supplied core and
    active orbitals, this includes the core-active cross-overlap blocks explicitly.
    """
    gto, fci, cistring = _imports()

    if previous.ncore != current.ncore:
        raise ValueError("Cannot overlap snapshots with different ncore.")
    if previous.ncas != current.ncas:
        raise ValueError("Cannot overlap snapshots with different ncas.")
    if tuple(previous.nelecas) != tuple(current.nelecas):
        raise ValueError("Cannot overlap snapshots with different active electron counts.")
    if previous.nroots != current.nroots:
        raise ValueError("Cannot track different numbers of roots.")

    ncore = previous.ncore
    ncas = previous.ncas
    nelecas = previous.nelecas
    norb = ncore + ncas
    neleca, nelecb = nelecas
    nelec_full = (ncore + neleca, ncore + nelecb)

    s_ao = gto.intor_cross(
        "int1e_ovlp_sph",
        previous.mol,
        current.mol,
    )

    C1 = previous.correlated_mo_coeff
    C2 = current.correlated_mo_coeff
    s_mo = C1.conj().T @ s_ao @ C2

    prev_full = [
        embed_active_ci_with_doubly_occupied_core(
            ci,
            ncore,
            ncas,
            nelecas,
            cistring=cistring,
        )
        for ci in previous.ci_roots
    ]

    curr_full = [
        embed_active_ci_with_doubly_occupied_core(
            ci,
            ncore,
            ncas,
            nelecas,
            cistring=cistring,
        )
        for ci in current.ci_roots
    ]

    O = np.zeros((previous.nroots, current.nroots), dtype=complex)

    for i, bra in enumerate(prev_full):
        for j, ket in enumerate(curr_full):
            O[i, j] = fci.addons.overlap(
                bra,
                ket,
                norb,
                nelec_full,
                s=s_mo,
            )

    return O


def correlated_orbital_cross_overlap(previous, current):
    """Return the core+active one-particle overlap matrix for diagnostics."""
    gto, _, _ = _imports()

    s_ao = gto.intor_cross(
        "int1e_ovlp_sph",
        previous.mol,
        current.mol,
    )

    return (
        previous.correlated_mo_coeff.conj().T
        @ s_ao
        @ current.correlated_mo_coeff
    )
