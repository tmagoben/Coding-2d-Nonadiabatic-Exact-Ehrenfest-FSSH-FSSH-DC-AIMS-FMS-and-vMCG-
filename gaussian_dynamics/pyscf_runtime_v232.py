"""Pinned PySCF 2.13.1 runtime, provenance, and real NAC evidence.

This module is intentionally fail closed.  A merely importable PySCF package is
not sufficient: both the installed distribution and imported module must report
exactly 2.13.1, the SA-CASSCF NAC API must import, and runtime provenance must be
fingerprinted before a real calculation is accepted.

The execution environment used for the release has a PID-namespace mismatch:
``os.getpid()`` is not represented by ``/proc/<pid>``, while ``/proc/self`` is
valid.  PySCF's Linux memory telemetry reads the former.  The guarded context
below substitutes the mathematically equivalent ``/proc/self/statm`` lookup only
when that exact ``FileNotFoundError`` is observed.  It does not change molecular
integrals, SCF/CASSCF tolerances, gradients, NACs, overlaps, or convergence gates.
"""

from contextlib import contextmanager
from dataclasses import asdict, dataclass
import base64
import hashlib
import importlib
from importlib import metadata
import importlib.util
import json
import math
import os
from pathlib import Path
import platform
import sys

import numpy as np

from .pyscf_nac_convention_v232 import (
    PYSCF_NAC_EMPIRICAL_MAPPING_V232,
    PYSCF_NAC_UPSTREAM_DOCUMENTATION_V232,
    PYSCF_REQUIRED_VERSION_V232,
)


RUNTIME_SCHEMA_V232 = "gnd-pyscf-runtime-evidence-v0.23.2"
_THREAD_ENVIRONMENT_KEYS_V232 = (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "BLIS_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
)


def _canonical_json_bytes_v232(value):
    return json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    ).encode("utf-8")


def _sha256_file_v232(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _distribution_inventory_sha256_v232(distribution_name):
    distribution = metadata.distribution(distribution_name)
    entries = []
    content_entries = []
    verified_file_count = 0
    verified_size_bytes = 0
    for item in distribution.files or ():
        file_hash = item.hash
        entries.append(
            {
                "path": str(item),
                "hash": None
                if file_hash is None
                else f"{file_hash.mode}:{file_hash.value}",
                "size": item.size,
            }
        )
        if file_hash is not None:
            path = distribution.locate_file(item)
            digest = hashlib.new(file_hash.mode)
            size_bytes = 0
            with Path(path).open("rb") as handle:
                for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                    digest.update(chunk)
                    size_bytes += len(chunk)
            encoded = base64.urlsafe_b64encode(digest.digest()).rstrip(b"=").decode()
            if encoded != file_hash.value:
                raise RuntimeError(
                    "installed PySCF file differs from distribution RECORD: "
                    f"{item}"
                )
            if item.size is not None and size_bytes != int(item.size):
                raise RuntimeError(
                    "installed PySCF file size differs from distribution RECORD: "
                    f"{item}"
                )
            content_entries.append(
                {
                    "path": str(item),
                    "hash": f"{file_hash.mode}:{encoded}",
                    "size": size_bytes,
                }
            )
            verified_file_count += 1
            verified_size_bytes += size_bytes
    entries.sort(key=lambda item: item["path"])
    content_entries.sort(key=lambda item: item["path"])
    return (
        hashlib.sha256(_canonical_json_bytes_v232(entries)).hexdigest(),
        len(entries),
        hashlib.sha256(_canonical_json_bytes_v232(content_entries)).hexdigest(),
        verified_file_count,
        verified_size_bytes,
    )


@dataclass(frozen=True)
class PySCFRuntimeProbeV232:
    installed: bool
    required_version: str
    distribution_version: str | None
    module_version: str | None
    exact_version: bool
    nac_api_available: bool
    usable: bool
    failure_reason: str | None

    def as_dict(self):
        return asdict(self)


def probe_pyscf_runtime_v232():
    """Probe without accepting absent, shadowed, or version-skewed runtimes."""
    try:
        installed = importlib.util.find_spec("pyscf") is not None
    except (ImportError, AttributeError, ValueError) as exc:
        return PySCFRuntimeProbeV232(
            installed=False,
            required_version=PYSCF_REQUIRED_VERSION_V232,
            distribution_version=None,
            module_version=None,
            exact_version=False,
            nac_api_available=False,
            usable=False,
            failure_reason=f"PySCF discovery failed: {type(exc).__name__}: {exc}",
        )

    if not installed:
        return PySCFRuntimeProbeV232(
            installed=False,
            required_version=PYSCF_REQUIRED_VERSION_V232,
            distribution_version=None,
            module_version=None,
            exact_version=False,
            nac_api_available=False,
            usable=False,
            failure_reason="PySCF is not installed.",
        )

    try:
        distribution_version = metadata.version("pyscf")
    except metadata.PackageNotFoundError:
        distribution_version = None

    try:
        pyscf = importlib.import_module("pyscf")
        module_version = str(getattr(pyscf, "__version__", "unknown"))
    except Exception as exc:  # an installed but unimportable runtime is unusable
        return PySCFRuntimeProbeV232(
            installed=True,
            required_version=PYSCF_REQUIRED_VERSION_V232,
            distribution_version=distribution_version,
            module_version=None,
            exact_version=False,
            nac_api_available=False,
            usable=False,
            failure_reason=f"PySCF import failed: {type(exc).__name__}: {exc}",
        )

    exact_version = bool(
        distribution_version == PYSCF_REQUIRED_VERSION_V232
        and module_version == PYSCF_REQUIRED_VERSION_V232
    )

    try:
        nac_module = importlib.import_module("pyscf.nac.sacasscf")
        nac_api_available = hasattr(nac_module, "NonAdiabaticCouplings")
    except (ImportError, AttributeError):
        nac_api_available = False

    if not exact_version:
        failure_reason = (
            "PySCF version mismatch: required "
            f"{PYSCF_REQUIRED_VERSION_V232}, distribution={distribution_version!r}, "
            f"module={module_version!r}."
        )
    elif not nac_api_available:
        failure_reason = "PySCF SA-CASSCF NAC API is unavailable."
    else:
        failure_reason = None

    return PySCFRuntimeProbeV232(
        installed=True,
        required_version=PYSCF_REQUIRED_VERSION_V232,
        distribution_version=distribution_version,
        module_version=module_version,
        exact_version=exact_version,
        nac_api_available=nac_api_available,
        usable=bool(exact_version and nac_api_available),
        failure_reason=failure_reason,
    )


def require_pyscf_runtime_v232():
    probe = probe_pyscf_runtime_v232()
    if not probe.installed:
        raise ImportError(
            "PySCF 2.13.1 is not installed; v0.23.2 real-runtime validation "
            "fails closed."
        )
    if not probe.exact_version:
        raise RuntimeError(probe.failure_reason)
    if not probe.nac_api_available:
        raise RuntimeError(probe.failure_reason)
    return probe


def _validated_memory_pair_v232(value):
    try:
        rss_mb, vms_mb = (float(value[0]), float(value[1]))
    except (TypeError, ValueError, IndexError) as exc:
        raise RuntimeError("PySCF memory telemetry returned an invalid value.") from exc
    if not math.isfinite(rss_mb) or not math.isfinite(vms_mb):
        raise RuntimeError("PySCF memory telemetry returned non-finite data.")
    if rss_mb < 0.0 or vms_mb < 0.0 or rss_mb > vms_mb:
        raise RuntimeError("PySCF memory telemetry returned inconsistent data.")
    return rss_mb, vms_mb


def _current_memory_via_proc_self_v232():
    if not sys.platform.startswith("linux"):
        raise RuntimeError("the /proc/self memory fallback is Linux-specific.")
    with Path("/proc/self/statm").open("r", encoding="ascii") as handle:
        fields = handle.readline().split()
    if len(fields) < 2:
        raise RuntimeError("/proc/self/statm did not contain vms and rss fields.")
    try:
        pages_vms, pages_rss = int(fields[0]), int(fields[1])
        page_size = int(os.sysconf("SC_PAGE_SIZE"))
    except (TypeError, ValueError, OSError) as exc:
        raise RuntimeError("could not parse /proc/self/statm memory telemetry.") from exc
    if pages_vms < 0 or pages_rss < 0 or page_size <= 0:
        raise RuntimeError("/proc/self/statm memory telemetry is invalid.")
    return _validated_memory_pair_v232(
        (pages_rss * page_size / 1e6, pages_vms * page_size / 1e6)
    )


@dataclass(frozen=True)
class PySCFRuntimeFingerprintV232:
    schema: str
    required_pyscf_version: str
    pyscf_distribution_version: str
    pyscf_module_version: str
    pyscf_module_sha256: str
    pyscf_nac_module_sha256: str
    pyscf_distribution_inventory_sha256: str
    pyscf_distribution_file_count: int
    pyscf_verified_content_sha256: str
    pyscf_verified_file_count: int
    pyscf_verified_size_bytes: int
    python_version: str
    python_implementation: str
    python_executable_sha256: str
    platform: str
    machine: str
    byteorder: str
    numpy_version: str
    scipy_version: str
    h5py_version: str
    thread_environment: dict
    memory_probe_mode: str
    environment_sha256: str

    def as_dict(self):
        return asdict(self)

    def validate(self):
        if self.required_pyscf_version != PYSCF_REQUIRED_VERSION_V232:
            raise ValueError("runtime fingerprint has the wrong required PySCF version.")
        if self.pyscf_distribution_version != PYSCF_REQUIRED_VERSION_V232:
            raise ValueError("runtime fingerprint has the wrong PySCF distribution version.")
        if self.pyscf_module_version != PYSCF_REQUIRED_VERSION_V232:
            raise ValueError("runtime fingerprint has the wrong imported PySCF version.")
        if self.memory_probe_mode not in {
            "pyscf_native_pid_statm",
            "proc_self_statm_pid_namespace_fallback",
            "proc_self_statm_requested",
        }:
            raise ValueError("runtime fingerprint has an unknown memory probe mode.")
        payload = self.as_dict()
        expected = payload.pop("environment_sha256")
        observed = hashlib.sha256(_canonical_json_bytes_v232(payload)).hexdigest()
        if observed != expected:
            raise ValueError("runtime environment fingerprint mismatch.")
        return self


def build_pyscf_runtime_fingerprint_v232(memory_probe_mode):
    probe = require_pyscf_runtime_v232()
    pyscf = importlib.import_module("pyscf")
    nac_module = importlib.import_module("pyscf.nac.sacasscf")
    (
        inventory_sha256,
        inventory_count,
        verified_content_sha256,
        verified_file_count,
        verified_size_bytes,
    ) = _distribution_inventory_sha256_v232("pyscf")
    payload = {
        "schema": "gnd-pyscf-runtime-fingerprint-v0.23.2",
        "required_pyscf_version": PYSCF_REQUIRED_VERSION_V232,
        "pyscf_distribution_version": str(probe.distribution_version),
        "pyscf_module_version": str(probe.module_version),
        "pyscf_module_sha256": _sha256_file_v232(pyscf.__file__),
        "pyscf_nac_module_sha256": _sha256_file_v232(nac_module.__file__),
        "pyscf_distribution_inventory_sha256": inventory_sha256,
        "pyscf_distribution_file_count": inventory_count,
        "pyscf_verified_content_sha256": verified_content_sha256,
        "pyscf_verified_file_count": verified_file_count,
        "pyscf_verified_size_bytes": verified_size_bytes,
        "python_version": platform.python_version(),
        "python_implementation": platform.python_implementation(),
        "python_executable_sha256": _sha256_file_v232(Path(sys.executable).resolve()),
        "platform": platform.platform(),
        "machine": platform.machine(),
        "byteorder": sys.byteorder,
        "numpy_version": metadata.version("numpy"),
        "scipy_version": metadata.version("scipy"),
        "h5py_version": metadata.version("h5py"),
        "thread_environment": {
            key: os.environ.get(key) for key in _THREAD_ENVIRONMENT_KEYS_V232
        },
        "memory_probe_mode": str(memory_probe_mode),
    }
    payload["environment_sha256"] = hashlib.sha256(
        _canonical_json_bytes_v232(payload)
    ).hexdigest()
    return PySCFRuntimeFingerprintV232(**payload).validate()


@dataclass(frozen=True)
class PySCFRuntimeContextV232:
    probe: PySCFRuntimeProbeV232
    fingerprint: PySCFRuntimeFingerprintV232
    memory_probe_mode: str
    memory_rss_mb_at_entry: float
    memory_vms_mb_at_entry: float


@contextmanager
def guarded_pyscf_runtime_v232(*, memory_probe_policy="auto"):
    """Yield an exact-version runtime with scoped memory telemetry.

    ``auto`` preserves native PySCF behavior unless the PID namespace requires the
    ``/proc/self`` fallback.  ``proc_self`` selects that equivalent Linux path
    deterministically for byte-reproducible release evidence.
    """
    if memory_probe_policy not in {"auto", "proc_self"}:
        raise ValueError("memory_probe_policy must be 'auto' or 'proc_self'.")
    probe = require_pyscf_runtime_v232()
    pyscf = importlib.import_module("pyscf")
    lib = pyscf.lib
    misc = importlib.import_module("pyscf.lib.misc")
    original_lib_current_memory = lib.current_memory
    original_misc_current_memory = misc.current_memory
    patched = False

    try:
        if memory_probe_policy == "proc_self":
            rss_mb, vms_mb = _current_memory_via_proc_self_v232()
            lib.current_memory = _current_memory_via_proc_self_v232
            misc.current_memory = _current_memory_via_proc_self_v232
            patched = True
            memory_probe_mode = "proc_self_statm_requested"
        else:
            try:
                rss_mb, vms_mb = _validated_memory_pair_v232(
                    original_lib_current_memory()
                )
                memory_probe_mode = "pyscf_native_pid_statm"
            except FileNotFoundError as exc:
                expected_path = f"/proc/{os.getpid()}/statm"
                if not sys.platform.startswith("linux") or exc.filename != expected_path:
                    raise
                rss_mb, vms_mb = _current_memory_via_proc_self_v232()
                lib.current_memory = _current_memory_via_proc_self_v232
                misc.current_memory = _current_memory_via_proc_self_v232
                patched = True
                memory_probe_mode = "proc_self_statm_pid_namespace_fallback"

        fingerprint = build_pyscf_runtime_fingerprint_v232(memory_probe_mode)
        yield PySCFRuntimeContextV232(
            probe=probe,
            fingerprint=fingerprint,
            memory_probe_mode=memory_probe_mode,
            memory_rss_mb_at_entry=rss_mb,
            memory_vms_mb_at_entry=vms_mb,
        )
    finally:
        if patched:
            lib.current_memory = original_lib_current_memory
            misc.current_memory = original_misc_current_memory


def _real_nested_v232(array):
    array = np.asarray(array, dtype=float)
    if array.ndim == 0:
        return float(array)
    return tuple(_real_nested_v232(item) for item in array)


def _complex_pairs_v232(array):
    array = np.asarray(array, dtype=complex)
    return _real_nested_v232(np.stack((array.real, array.imag), axis=-1))


def _phase_aligned_overlap_v232(reference, displaced):
    from .pyscf_wavefunction_overlap import casscf_state_overlap_matrix

    overlap = np.asarray(
        casscf_state_overlap_matrix(reference, displaced), dtype=complex
    )
    diagonal = np.diag(overlap)
    if np.any(np.abs(diagonal) < 0.8):
        raise RuntimeError("displaced SA-CASSCF roots cannot be phase aligned safely.")
    phases = np.conj(diagonal) / np.abs(diagonal)
    return overlap * phases[None, :]


def _h3p_geometry_v232(coords_bohr):
    from .molecular_backend import MolecularGeometry

    return MolecularGeometry(("H", "H", "H"), np.asarray(coords_bohr, dtype=float))


def _h3p_config_v232(*, use_etfs, compute_scaled_nac=False):
    from .pyscf_backend_v05 import PySCFSACASSCFConfig

    return PySCFSACASSCFConfig(
        basis="sto-3g",
        ncas=3,
        nelecas=2,
        nstates=3,
        weights=(1.0 / 3.0,) * 3,
        charge=1,
        spin=0,
        symmetry=False,
        scf_reference="RHF",
        scf_conv_tol=1e-12,
        scf_max_cycle=100,
        mc_conv_tol=1e-11,
        mc_conv_tol_grad=1e-7,
        mc_max_cycle_macro=50,
        use_etfs=bool(use_etfs),
        compute_scaled_nac=bool(compute_scaled_nac),
        warm_start_mo=False,
        verbose=0,
        max_memory_mb=1000,
    )


def _production_h3p_point_v232(coords_bohr, *, use_etfs, compute_scaled_nac=False):
    from .pyscf_tracked_backend_v06 import PySCFTrackedSACASSCFBackend

    backend = PySCFTrackedSACASSCFBackend(
        _h3p_config_v232(
            use_etfs=use_etfs,
            compute_scaled_nac=compute_scaled_nac,
        )
    )
    return backend.evaluate_raw_with_snapshot(_h3p_geometry_v232(coords_bohr))


def _direct_h3p_solver_v232(coords_bohr):
    """Run a converged SA-CASSCF solver retained for tuple-level diagnostics."""
    from pyscf import gto, mcscf, scf

    from .pyscf_wavefunction_overlap import CASSCFWavefunctionSnapshot

    coords_bohr = np.asarray(coords_bohr, dtype=float)
    mol = gto.M(
        atom=[
            (symbol, tuple(float(x) for x in row))
            for symbol, row in zip(("H", "H", "H"), coords_bohr)
        ],
        basis="sto-3g",
        charge=1,
        spin=0,
        unit="Bohr",
        symmetry=False,
        verbose=0,
        max_memory=1000,
    )
    mf = scf.RHF(mol)
    mf.conv_tol = 1e-12
    mf.max_cycle = 100
    mf.kernel()
    if not bool(mf.converged):
        raise RuntimeError("direct H3+ PySCF RHF did not converge.")

    mc = mcscf.CASSCF(mf, 3, 2).state_average_((1.0 / 3.0,) * 3)
    mc.conv_tol = 1e-11
    mc.conv_tol_grad = 1e-7
    mc.max_cycle_macro = 50
    mc.kernel()
    if not bool(mc.converged):
        raise RuntimeError("direct H3+ PySCF SA-CASSCF did not converge.")

    snapshot = CASSCFWavefunctionSnapshot(
        mol=mol,
        mo_coeff=np.asarray(mc.mo_coeff).copy(),
        ci_roots=tuple(np.asarray(root).copy() for root in mc.ci),
        ncore=int(mc.ncore),
        ncas=int(mc.ncas),
        nelecas=tuple(int(x) for x in mc.nelecas),
        metadata={
            "runtime": "PySCF",
            "runtime_version": PYSCF_REQUIRED_VERSION_V232,
            "calculation": "H3+ SA-CASSCF tuple-orientation diagnostic",
        },
    )
    return mc, snapshot


@dataclass(frozen=True)
class PySCFRuntimeEvidenceV232:
    schema: str
    runtime: PySCFRuntimeFingerprintV232
    smoke: dict
    nac_mapping: dict
    checks: dict

    @property
    def passed(self):
        return bool(self.checks and all(value is True for value in self.checks.values()))

    def as_dict(self):
        payload = asdict(self)
        payload["passed"] = self.passed
        payload["evidence_sha256"] = self.fingerprint()
        return payload

    def fingerprint(self):
        payload = asdict(self)
        return hashlib.sha256(_canonical_json_bytes_v232(payload)).hexdigest()

    def validate(self):
        self.runtime.validate()
        if self.schema != RUNTIME_SCHEMA_V232:
            raise ValueError("unknown PySCF runtime evidence schema.")
        if not self.passed:
            failed = [name for name, value in self.checks.items() if value is not True]
            raise RuntimeError("PySCF runtime evidence failed: " + ", ".join(failed))
        return self


def run_pyscf_runtime_evidence_v232(
    step_sizes=(1e-2, 1e-3, 1e-4),
    *,
    memory_probe_policy="auto",
):
    """Run genuine energy/gradient/NAC/overlap and NAC-sign certification.

    The molecule is asymmetric triangular H3+ at SA-CASSCF(2e,3o)/STO-3G with
    three equally weighted roots.  Its 0/2 derivative coupling is nonzero, which
    makes the central-difference sign test discriminating rather than ceremonial.
    """
    steps = tuple(float(value) for value in step_sizes)
    if len(steps) < 3 or any(value <= 0.0 for value in steps):
        raise ValueError("at least three positive central-difference steps are required.")
    if any(steps[index + 1] >= steps[index] for index in range(len(steps) - 1)):
        raise ValueError("central-difference steps must be strictly decreasing.")

    center = np.array(
        [[0.0, 0.0, 0.0], [1.5, 0.0, 0.0], [0.2, 1.3, 0.0]],
        dtype=float,
    )
    atom_index = 2
    cartesian_index = 1
    state_i, state_j = 0, 2

    with guarded_pyscf_runtime_v232(
        memory_probe_policy=memory_probe_policy
    ) as runtime:
        point_center, snapshot_center = _production_h3p_point_v232(
            center, use_etfs=False, compute_scaled_nac=True
        )
        point_etf, _ = _production_h3p_point_v232(
            center, use_etfs=True, compute_scaled_nac=False
        )

        right_coords = center.copy()
        right_coords[atom_index, cartesian_index] += steps[0]
        point_right, snapshot_right = _production_h3p_point_v232(
            right_coords, use_etfs=False, compute_scaled_nac=False
        )

        direct_mc, direct_center = _direct_h3p_solver_v232(center)
        direct_nac = direct_mc.nac_method()
        raw_no_etf_ij = np.asarray(
            direct_nac.kernel(
                state=(state_i, state_j), use_etfs=False, mult_ediff=False
            ),
            dtype=float,
        )
        raw_no_etf_ji = np.asarray(
            direct_nac.kernel(
                state=(state_j, state_i), use_etfs=False, mult_ediff=False
            ),
            dtype=float,
        )
        raw_etf_ij = np.asarray(
            direct_nac.kernel(
                state=(state_i, state_j), use_etfs=True, mult_ediff=False
            ),
            dtype=float,
        )
        raw_etf_ji = np.asarray(
            direct_nac.kernel(
                state=(state_j, state_i), use_etfs=True, mult_ediff=False
            ),
            dtype=float,
        )
        raw_scaled_ij = np.asarray(
            direct_nac.kernel(
                state=(state_i, state_j), use_etfs=False, mult_ediff=True
            ),
            dtype=float,
        )
        raw_scaled_ji = np.asarray(
            direct_nac.kernel(
                state=(state_j, state_i), use_etfs=False, mult_ediff=True
            ),
            dtype=float,
        )

        production_fd = []
        direct_fd = []
        production_errors = []
        direct_errors = []
        displaced_snapshots = {}
        for step in steps:
            if step == steps[0]:
                plus_snapshot = snapshot_right
            else:
                plus_coords = center.copy()
                plus_coords[atom_index, cartesian_index] += step
                _, plus_snapshot = _direct_h3p_solver_v232(plus_coords)

            minus_coords = center.copy()
            minus_coords[atom_index, cartesian_index] -= step
            _, minus_snapshot = _direct_h3p_solver_v232(minus_coords)
            displaced_snapshots[step] = (plus_snapshot, minus_snapshot)

            plus_production = _phase_aligned_overlap_v232(
                snapshot_center, plus_snapshot
            )
            minus_production = _phase_aligned_overlap_v232(
                snapshot_center, minus_snapshot
            )
            fd_production = (plus_production - minus_production) / (2.0 * step)
            production_fd.append(fd_production)
            target = point_center.nac_cart[:, :, atom_index, cartesian_index]
            production_errors.append(float(np.max(np.abs(fd_production - target))))

            plus_direct = _phase_aligned_overlap_v232(direct_center, plus_snapshot)
            minus_direct = _phase_aligned_overlap_v232(direct_center, minus_snapshot)
            fd_direct = (plus_direct - minus_direct) / (2.0 * step)
            direct_fd.append(fd_direct)
            direct_errors.append(
                float(
                    abs(
                        fd_direct[state_i, state_j]
                        - raw_no_etf_ij[atom_index, cartesian_index]
                    )
                )
            )

        from .pyscf_wavefunction_overlap import casscf_state_overlap_matrix

        overlap_lr = np.asarray(
            casscf_state_overlap_matrix(snapshot_center, snapshot_right),
            dtype=complex,
        )
        overlap_rl = np.asarray(
            casscf_state_overlap_matrix(snapshot_right, snapshot_center),
            dtype=complex,
        )
        overlap_ll = np.asarray(
            casscf_state_overlap_matrix(snapshot_center, snapshot_center),
            dtype=complex,
        )
        overlap_rr = np.asarray(
            casscf_state_overlap_matrix(snapshot_right, snapshot_right),
            dtype=complex,
        )
        singular_values = np.linalg.svd(overlap_lr, compute_uv=False)
        cross_isometry_defect = float(
            np.linalg.norm(overlap_lr.conj().T @ overlap_lr - np.eye(3), ord="fro")
        )
        reciprocity_residual = float(
            np.linalg.norm(overlap_lr - overlap_rl.conj().T, ord="fro")
        )
        self_overlap_residual = float(
            max(
                np.linalg.norm(overlap_ll - np.eye(3), ord="fro"),
                np.linalg.norm(overlap_rr - np.eye(3), ord="fro"),
            )
        )
        contraction_excess = float(max(0.0, np.max(singular_values) - 1.0))

        nac_stack = np.stack((point_center.nac_cart, point_right.nac_cart))
        gradient_stack = np.stack(
            (point_center.gradients_cart, point_right.gradients_cart)
        )
        energy_stack = np.stack((point_center.energies, point_right.energies))
        nac_antisymmetry_residual = float(
            np.max(np.abs(nac_stack + np.swapaxes(nac_stack, 1, 2)))
        )
        gradient_translation_residual = float(
            np.max(np.abs(np.sum(gradient_stack, axis=2)))
        )
        no_etf_translation_residual = float(
            np.max(np.abs(np.sum(raw_no_etf_ij, axis=0)))
        )
        etf_translation_residual = float(
            np.max(np.abs(np.sum(raw_etf_ij, axis=0)))
        )
        tuple_antisymmetry_no_etf = float(
            np.max(np.abs(raw_no_etf_ij + raw_no_etf_ji))
        )
        tuple_antisymmetry_etf = float(np.max(np.abs(raw_etf_ij + raw_etf_ji)))
        scaled_tuple_symmetry = float(np.max(np.abs(raw_scaled_ij - raw_scaled_ji)))
        direct_energy_gap = float(
            np.asarray(direct_mc.e_states)[state_j]
            - np.asarray(direct_mc.e_states)[state_i]
        )
        scaled_energy_relation = float(
            np.max(np.abs(raw_scaled_ij - direct_energy_gap * raw_no_etf_ij))
        )
        if point_center.scaled_nac_cart is None:
            raise RuntimeError("production scaled-NAC diagnostic was not produced.")
        production_scaled_symmetry = float(
            np.max(
                np.abs(
                    point_center.scaled_nac_cart
                    - np.swapaxes(point_center.scaled_nac_cart, 0, 1)
                )
            )
        )
        production_scaled_energy_relation = 0.0
        for left_state in range(3):
            for right_state in range(left_state + 1, 3):
                gap = (
                    point_center.energies[right_state]
                    - point_center.energies[left_state]
                )
                production_scaled_energy_relation = max(
                    production_scaled_energy_relation,
                    float(
                        np.max(
                            np.abs(
                                point_center.scaled_nac_cart[
                                    left_state, right_state
                                ]
                                - gap
                                * point_center.nac_cart[left_state, right_state]
                            )
                        )
                    ),
                )

        smoke = {
            "system": "H3+ asymmetric triangle",
            "symbols": ["H", "H", "H"],
            "coordinates_bohr": _real_nested_v232(
                np.stack((center, right_coords))
            ),
            "method": "SA-CASSCF(2e,3o)/STO-3G; three equal-weight singlet roots",
            "energies_hartree": _real_nested_v232(energy_stack),
            "gradients_hartree_per_bohr": _real_nested_v232(gradient_stack),
            "nac_inverse_bohr": _real_nested_v232(nac_stack),
            "overlap_left_right": _complex_pairs_v232(overlap_lr),
            "overlap_right_left": _complex_pairs_v232(overlap_rl),
            "overlap_left_left": _complex_pairs_v232(overlap_ll),
            "overlap_right_right": _complex_pairs_v232(overlap_rr),
            "singular_values": _real_nested_v232(singular_values),
            "metrics": {
                "nac_antisymmetry_residual": nac_antisymmetry_residual,
                "gradient_translation_residual": gradient_translation_residual,
                "nac_norm": float(np.linalg.norm(nac_stack)),
                "gradient_norm": float(np.linalg.norm(gradient_stack)),
                "minimum_energy_gap_hartree": float(
                    np.min(np.diff(energy_stack, axis=1))
                ),
                "energy_response_hartree": float(
                    np.max(np.abs(energy_stack[1] - energy_stack[0]))
                ),
                "reciprocity_residual": reciprocity_residual,
                "self_overlap_residual": self_overlap_residual,
                "contraction_excess": contraction_excess,
                "minimum_singular_value": float(np.min(singular_values)),
                "cross_overlap_isometry_defect": cross_isometry_defect,
            },
        }

        nac_mapping = {
            "upstream_documentation": PYSCF_NAC_UPSTREAM_DOCUMENTATION_V232,
            "empirical_mapping": PYSCF_NAC_EMPIRICAL_MAPPING_V232,
            "state_pair": [state_i, state_j],
            "displaced_coordinate": {
                "atom_index": atom_index,
                "cartesian_index": cartesian_index,
                "label": "H[2] y",
            },
            "central_difference_steps_bohr": list(steps),
            "production_target_matrix_inverse_bohr": _real_nested_v232(
                point_center.nac_cart[:, :, atom_index, cartesian_index]
            ),
            "production_central_difference_matrices_inverse_bohr": tuple(
                _complex_pairs_v232(value) for value in production_fd
            ),
            "direct_central_difference_matrices_inverse_bohr": tuple(
                _complex_pairs_v232(value) for value in direct_fd
            ),
            "raw_state_i_j_no_etf_inverse_bohr": _real_nested_v232(raw_no_etf_ij),
            "raw_state_j_i_no_etf_inverse_bohr": _real_nested_v232(raw_no_etf_ji),
            "raw_state_i_j_etf_inverse_bohr": _real_nested_v232(raw_etf_ij),
            "raw_state_j_i_etf_inverse_bohr": _real_nested_v232(raw_etf_ji),
            "raw_scaled_state_i_j_hartree_per_bohr": _real_nested_v232(raw_scaled_ij),
            "raw_scaled_state_j_i_hartree_per_bohr": _real_nested_v232(raw_scaled_ji),
            "metrics": {
                "production_central_difference_max_errors": production_errors,
                "direct_state_i_j_central_difference_errors": direct_errors,
                "production_dij_selected": float(
                    point_center.nac_cart[
                        state_i, state_j, atom_index, cartesian_index
                    ]
                ),
                "finest_production_fd_selected_real": float(
                    np.real(production_fd[-1][state_i, state_j])
                ),
                "raw_state_i_j_selected": float(
                    raw_no_etf_ij[atom_index, cartesian_index]
                ),
                "raw_state_j_i_selected": float(
                    raw_no_etf_ji[atom_index, cartesian_index]
                ),
                "tuple_antisymmetry_no_etf": tuple_antisymmetry_no_etf,
                "tuple_antisymmetry_etf": tuple_antisymmetry_etf,
                "scaled_tuple_symmetry": scaled_tuple_symmetry,
                "scaled_energy_relation_residual": scaled_energy_relation,
                "production_scaled_tuple_symmetry": production_scaled_symmetry,
                "production_scaled_energy_relation_residual": (
                    production_scaled_energy_relation
                ),
                "no_etf_translation_residual": no_etf_translation_residual,
                "etf_translation_residual": etf_translation_residual,
                "etf_vs_no_etf_norm": float(
                    np.linalg.norm(raw_etf_ij - raw_no_etf_ij)
                ),
            },
        }

        production_converges_quadratically = all(
            production_errors[index + 1] < 0.05 * production_errors[index]
            for index in range(len(production_errors) - 1)
        )
        direct_converges_quadratically = all(
            direct_errors[index + 1] < 0.05 * direct_errors[index]
            for index in range(len(direct_errors) - 1)
        )
        all_numeric = all(
            np.all(np.isfinite(value))
            for value in (energy_stack, gradient_stack, nac_stack, overlap_lr)
        )
        checks = {
            "runtime_exactly_pyscf_2_13_1": runtime.probe.usable is True,
            "runtime_fingerprint_valid": (
                runtime.fingerprint.validate() is runtime.fingerprint
            ),
            "sandbox_memory_probe_is_valid": (
                runtime.memory_rss_mb_at_entry >= 0.0
                and runtime.memory_vms_mb_at_entry >= runtime.memory_rss_mb_at_entry
            ),
            "two_geometry_values_are_finite": bool(all_numeric),
            "two_geometry_energy_response_is_nonzero": (
                smoke["metrics"]["energy_response_hartree"] > 1e-8
            ),
            "analytic_gradients_are_nontrivial": (
                smoke["metrics"]["gradient_norm"] > 1e-6
            ),
            "analytic_gradient_translation_invariance": (
                gradient_translation_residual < 1e-9
            ),
            "nac_is_nontrivial": smoke["metrics"]["nac_norm"] > 1e-6,
            "nac_is_antisymmetric": nac_antisymmetry_residual < 1e-12,
            "energy_roots_are_separated": (
                smoke["metrics"]["minimum_energy_gap_hartree"] > 1e-4
            ),
            "self_overlap_is_identity": self_overlap_residual < 1e-10,
            "cross_overlap_is_reciprocal": reciprocity_residual < 1e-10,
            "cross_overlap_is_a_physical_contraction": contraction_excess < 1e-10,
            "represented_manifold_overlap_is_retained": (
                smoke["metrics"]["minimum_singular_value"] > 0.99
            ),
            "cross_overlap_is_not_forced_to_exact_isometry": (
                cross_isometry_defect > 1e-8
            ),
            "production_state_tuple_mapping_has_correct_sign": (
                production_errors[-1] < 1e-6
            ),
            "production_mapping_shows_second_order_convergence": (
                production_converges_quadratically
            ),
            "raw_state_i_j_matches_overlap_derivative": direct_errors[-1] < 1e-6,
            "raw_state_i_j_mapping_shows_second_order_convergence": (
                direct_converges_quadratically
            ),
            "opposite_state_tuples_are_antisymmetric_no_etf": (
                tuple_antisymmetry_no_etf < 1e-10
            ),
            "opposite_state_tuples_are_antisymmetric_etf": (
                tuple_antisymmetry_etf < 1e-10
            ),
            "scaled_nac_is_symmetric_under_tuple_swap": (
                scaled_tuple_symmetry < 1e-10
            ),
            "scaled_nac_obeys_energy_difference_relation": (
                scaled_energy_relation < 1e-10
            ),
            "production_scaled_nac_is_symmetric": (
                production_scaled_symmetry < 1e-10
            ),
            "production_scaled_nac_obeys_energy_difference_relation": (
                production_scaled_energy_relation < 1e-10
            ),
            "full_overlap_derivative_uses_no_etf": (
                no_etf_translation_residual > 1e-3
            ),
            "etf_removes_translation_component": etf_translation_residual < 1e-10,
            "etf_and_full_overlap_derivatives_are_distinct": (
                nac_mapping["metrics"]["etf_vs_no_etf_norm"] > 1e-3
            ),
        }

        return PySCFRuntimeEvidenceV232(
            schema=RUNTIME_SCHEMA_V232,
            runtime=runtime.fingerprint,
            smoke=smoke,
            nac_mapping=nac_mapping,
            checks={name: bool(value) for name, value in checks.items()},
        ).validate()
