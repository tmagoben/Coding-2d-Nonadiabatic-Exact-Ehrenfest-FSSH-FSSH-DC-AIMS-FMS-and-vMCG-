"""Separate byte-reproducible release identity from scientific compatibility."""

from dataclasses import asdict, dataclass, is_dataclass
import hashlib
import json
import math


RUNTIME_PROFILE_SCHEMA_V233 = "gnd-runtime-compatibility-profile-v0.23.3"
RUNTIME_PROFILE_KINDS_V233 = ("release_locked", "scientifically_compatible")


def _version_tuple_v233(value):
    pieces = str(value).split(".")
    if not pieces or any(not piece.isdigit() for piece in pieces):
        raise ValueError(f"version {value!r} is not a dotted numeric version.")
    return tuple(int(piece) for piece in pieces)


def _canonical_bytes_v233(value):
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")


@dataclass(frozen=True)
class RuntimeCompatibilityProfileV233:
    schema: str
    kind: str
    required_pyscf_version: str
    python_minimum: str
    python_maximum_exclusive: str
    numpy_minimum: str
    numpy_maximum_exclusive: str
    scipy_minimum: str
    scipy_maximum_exclusive: str
    h5py_minimum: str
    h5py_maximum_exclusive: str
    exact_identity: dict | None = None

    def validate(self):
        if self.schema != RUNTIME_PROFILE_SCHEMA_V233:
            raise ValueError("runtime profile schema mismatch.")
        if self.kind not in RUNTIME_PROFILE_KINDS_V233:
            raise ValueError("runtime profile kind is unsupported.")
        for name in (
            "required_pyscf_version",
            "python_minimum",
            "python_maximum_exclusive",
            "numpy_minimum",
            "numpy_maximum_exclusive",
            "scipy_minimum",
            "scipy_maximum_exclusive",
            "h5py_minimum",
            "h5py_maximum_exclusive",
        ):
            _version_tuple_v233(getattr(self, name))
        if self.kind == "release_locked":
            if not isinstance(self.exact_identity, dict) or not self.exact_identity:
                raise ValueError("release-locked profile requires exact identity.")
        elif self.exact_identity is not None:
            raise ValueError(
                "scientific-compatibility profile cannot contain exact identity."
            )
        return self

    def as_dict(self):
        return asdict(self)

    def fingerprint(self):
        self.validate()
        return hashlib.sha256(_canonical_bytes_v233(self.as_dict())).hexdigest()


_RELEASE_EXACT_IDENTITY_V233 = {
    "python_version": "3.12.13",
    "python_implementation": "CPython",
    "machine": "x86_64",
    "byteorder": "little",
    "numpy_version": "2.5.2",
    "scipy_version": "1.18.0",
    "h5py_version": "3.16.0",
    "pyscf_distribution_version": "2.13.1",
    "pyscf_module_version": "2.13.1",
    "pyscf_module_sha256": (
        "9c06579656d120a595ce5b84a5066aa78d6eec547a7c6a6388799f909a461f92"
    ),
    "pyscf_nac_module_sha256": (
        "48737e950a469a56b485af58a571ac9f588c16605a39b60f78abb1bb96aa7c91"
    ),
    "pyscf_verified_content_sha256": (
        "962f2a5ff9e071bda6785bf89b169938a4e4f1adf81347b942ffa32c8b0c5347"
    ),
    "pyscf_verified_file_count": 1193,
    "pyscf_verified_size_bytes": 168868242,
    "python_executable_sha256": (
        "021044895e95be79dc2f110367607e684119afbc8ce75f6f0eec94844e0acec7"
    ),
    "thread_environment": {
        "BLIS_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
        "NUMEXPR_NUM_THREADS": "1",
        "OMP_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
        "VECLIB_MAXIMUM_THREADS": "1",
    },
}


def release_locked_runtime_profile_v233():
    return RuntimeCompatibilityProfileV233(
        schema=RUNTIME_PROFILE_SCHEMA_V233,
        kind="release_locked",
        required_pyscf_version="2.13.1",
        python_minimum="3.12.13",
        python_maximum_exclusive="3.12.14",
        numpy_minimum="2.5.2",
        numpy_maximum_exclusive="2.5.3",
        scipy_minimum="1.18.0",
        scipy_maximum_exclusive="1.18.1",
        h5py_minimum="3.16.0",
        h5py_maximum_exclusive="3.16.1",
        exact_identity=dict(_RELEASE_EXACT_IDENTITY_V233),
    ).validate()


def scientifically_compatible_runtime_profile_v233():
    return RuntimeCompatibilityProfileV233(
        schema=RUNTIME_PROFILE_SCHEMA_V233,
        kind="scientifically_compatible",
        required_pyscf_version="2.13.1",
        python_minimum="3.10.0",
        python_maximum_exclusive="3.14.0",
        numpy_minimum="1.24.0",
        numpy_maximum_exclusive="3.0.0",
        scipy_minimum="1.10.0",
        scipy_maximum_exclusive="2.0.0",
        h5py_minimum="3.8.0",
        h5py_maximum_exclusive="4.0.0",
        exact_identity=None,
    ).validate()


@dataclass(frozen=True)
class RuntimeCompatibilityReportV233:
    profile_kind: str
    profile_fingerprint: str
    compatible: bool
    checks: dict
    mismatches: tuple

    def as_dict(self):
        return {
            "profile_kind": self.profile_kind,
            "profile_fingerprint": self.profile_fingerprint,
            "compatible": bool(self.compatible),
            "checks": dict(self.checks),
            "mismatches": list(self.mismatches),
        }


def _runtime_mapping_v233(runtime):
    if is_dataclass(runtime):
        return asdict(runtime)
    if isinstance(runtime, dict):
        return dict(runtime)
    if hasattr(runtime, "as_dict"):
        payload = runtime.as_dict()
        if isinstance(payload, dict):
            return dict(payload)
    raise TypeError("runtime evidence must be a dataclass or mapping.")


def _in_range_v233(value, minimum, maximum_exclusive):
    observed = _version_tuple_v233(value)
    return bool(
        _version_tuple_v233(minimum)
        <= observed
        < _version_tuple_v233(maximum_exclusive)
    )


def assess_runtime_compatibility_v233(runtime, profile):
    """Assess exact release identity or the broader declared scientific range."""
    if type(profile) is not RuntimeCompatibilityProfileV233:
        raise TypeError("profile must be RuntimeCompatibilityProfileV233.")
    profile = profile.validate()
    runtime = _runtime_mapping_v233(runtime)
    checks = {
        "pyscf_distribution_exact": runtime.get("pyscf_distribution_version")
        == profile.required_pyscf_version,
        "pyscf_module_exact": runtime.get("pyscf_module_version")
        == profile.required_pyscf_version,
        "python_in_range": _in_range_v233(
            runtime.get("python_version", "0"),
            profile.python_minimum,
            profile.python_maximum_exclusive,
        ),
        "numpy_in_range": _in_range_v233(
            runtime.get("numpy_version", "0"),
            profile.numpy_minimum,
            profile.numpy_maximum_exclusive,
        ),
        "scipy_in_range": _in_range_v233(
            runtime.get("scipy_version", "0"),
            profile.scipy_minimum,
            profile.scipy_maximum_exclusive,
        ),
        "h5py_in_range": _in_range_v233(
            runtime.get("h5py_version", "0"),
            profile.h5py_minimum,
            profile.h5py_maximum_exclusive,
        ),
        "supported_byteorder": runtime.get("byteorder") == "little",
        "supported_machine": runtime.get("machine") in {"x86_64", "aarch64"},
    }
    if profile.kind == "release_locked":
        for name, expected in profile.exact_identity.items():
            checks[f"exact::{name}"] = runtime.get(name) == expected
        platform_value = str(runtime.get("platform", ""))
        checks["exact::linux_platform"] = platform_value.startswith("Linux-")
        checks["exact::memory_probe"] = runtime.get("memory_probe_mode") in {
            "proc_self_statm_requested",
            "native_pid_statm",
            "proc_self_statm_pid_namespace_fallback",
        }
    checks = {name: bool(value) for name, value in checks.items()}
    mismatches = tuple(name for name, value in checks.items() if not value)
    return RuntimeCompatibilityReportV233(
        profile_kind=profile.kind,
        profile_fingerprint=profile.fingerprint(),
        compatible=not mismatches,
        checks=checks,
        mismatches=mismatches,
    )
