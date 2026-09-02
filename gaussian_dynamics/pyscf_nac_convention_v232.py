"""Empirically certified PySCF 2.13.1 SA-CASSCF NAC mapping.

The upstream PySCF docstring describes ``state=(ket, bra)`` as returning
``<bra|d ket/dR>``.  For PySCF 2.13.1, phase-aligned central differences of
independently calculated many-electron SA-CASSCF overlaps instead show that the
array returned for ``state=(a, b)`` has the orientation ``<a|d b/dR>``.

v0.23.2 therefore freezes the *empirically validated* production mapping below.
This constant is deliberately distinct from the upstream documentation string;
silently treating those two statements as equivalent caused the pre-v0.23.2
production paths to populate the internal derivative-coupling matrix with the
wrong sign.
"""

from importlib import metadata


PYSCF_REQUIRED_VERSION_V232 = "2.13.1"

PYSCF_NAC_UPSTREAM_DOCUMENTATION_V232 = (
    "state=(ket,bra) is documented as returning <bra|d ket/dR>"
)

PYSCF_NAC_EMPIRICAL_MAPPING_V232 = (
    "PySCF 2.13.1 state=(i,j) populates internal d[i,j]=<i|d j/dR>; "
    "certified by phase-aligned many-electron overlap central differences "
    "with use_etfs=False"
)


def require_exact_pyscf_version_v232(
    pyscf_module,
    *,
    distribution_version=None,
):
    """Reject unpinned or shadowed PySCF production runtimes.

    Unit tests may inject ``distribution_version`` explicitly.  Production callers
    leave it unset so the installed distribution record is checked independently
    from the imported module.
    """
    observed = str(getattr(pyscf_module, "__version__", "unknown"))
    if distribution_version is None:
        try:
            distribution_version = metadata.version("pyscf")
        except metadata.PackageNotFoundError as exc:
            raise RuntimeError(
                "PySCF is importable but its installed distribution provenance "
                "is absent."
            ) from exc
    distribution_version = str(distribution_version)
    if (
        observed != PYSCF_REQUIRED_VERSION_V232
        or distribution_version != PYSCF_REQUIRED_VERSION_V232
    ):
        raise RuntimeError(
            "Gaussian Nonadiabatic Dynamics v0.23.2 requires exactly PySCF "
            f"{PYSCF_REQUIRED_VERSION_V232}; imported module={observed}, "
            f"distribution={distribution_version}."
        )
    return observed


def pyscf_state_tuple_for_internal_dij_v232(i, j):
    """Return the PySCF 2.13.1 state tuple for internal ``d[i,j]``.

    ``i == j`` is excluded because diagonal derivative couplings are fixed to
    zero by the real adiabatic-state gauge used by these production backends.
    """
    i = int(i)
    j = int(j)
    if i < 0 or j < 0:
        raise ValueError("state indices must be nonnegative.")
    if i == j:
        raise ValueError("off-diagonal state indices must differ.")
    return i, j
