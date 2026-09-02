"""Derivative-coupling convention identity and legacy-data quarantine for v0.23.3."""

from dataclasses import asdict, dataclass
import hashlib
import json

from .pyscf_nac_convention_v232 import (
    PYSCF_NAC_EMPIRICAL_MAPPING_V232,
    PYSCF_REQUIRED_VERSION_V232,
)


NAC_CONVENTION_SCHEMA_V233 = "gnd-derivative-coupling-convention-v0.23.3"
INTERNAL_NAC_DEFINITION_V233 = "d[i,j]=<Phi_i|d Phi_j/dR>"
PYSCF_NAC_MAPPING_ID_V233 = "pyscf-2.13.1-state-i-j-no-etf-unscaled-v1"
PYSCF_ETF_NAC_MAPPING_ID_V233 = (
    "pyscf-2.13.1-state-i-j-etf-unscaled-diagnostic-v1"
)
ANALYTIC_NAC_MAPPING_ID_V233 = "analytic-provider-native-internal-dij-v1"

LEGACY_NAC_DISPOSITIONS_V233 = (
    "not_pyscf_derived",
    "verified_corrected_v232",
    "requires_sign_correction",
    "unknown",
)


def _canonical_json_bytes_v233(value):
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")


@dataclass(frozen=True)
class DerivativeCouplingConventionV233:
    schema: str
    internal_definition: str
    source_backend: str
    source_backend_version: str
    source_mapping_id: str
    source_mapping_description: str
    use_etfs: bool
    mult_ediff: bool
    coordinate_unit: str = "bohr"

    def validate(self):
        if self.schema != NAC_CONVENTION_SCHEMA_V233:
            raise ValueError("derivative-coupling convention schema mismatch.")
        if self.internal_definition != INTERNAL_NAC_DEFINITION_V233:
            raise ValueError("internal derivative-coupling definition mismatch.")
        for name in (
            "source_backend",
            "source_backend_version",
            "source_mapping_id",
            "source_mapping_description",
        ):
            if not str(getattr(self, name)).strip():
                raise ValueError(f"{name} cannot be empty.")
        if type(self.use_etfs) is not bool or type(self.mult_ediff) is not bool:
            raise TypeError("ETF and energy-scaling flags must be native Booleans.")
        if self.coordinate_unit != "bohr":
            raise ValueError("v0.23.3 derivative couplings require bohr coordinates.")
        if self.source_backend == "PySCF":
            if self.source_backend_version != PYSCF_REQUIRED_VERSION_V232:
                raise ValueError("PySCF NAC convention requires exactly version 2.13.1.")
            expected_mapping = (
                PYSCF_ETF_NAC_MAPPING_ID_V233
                if self.use_etfs
                else PYSCF_NAC_MAPPING_ID_V233
            )
            if self.source_mapping_id != expected_mapping:
                raise ValueError("PySCF NAC mapping identity and ETF policy disagree.")
            if self.mult_ediff:
                raise ValueError(
                    "dynamics derivative couplings require mult_ediff=False."
                )
        return self

    def as_dict(self):
        return asdict(self)

    def fingerprint(self):
        self.validate()
        return hashlib.sha256(_canonical_json_bytes_v233(self.as_dict())).hexdigest()


def corrected_pyscf_nac_convention_v233(*, use_etfs=False):
    use_etfs = bool(use_etfs)
    return DerivativeCouplingConventionV233(
        schema=NAC_CONVENTION_SCHEMA_V233,
        internal_definition=INTERNAL_NAC_DEFINITION_V233,
        source_backend="PySCF",
        source_backend_version=PYSCF_REQUIRED_VERSION_V232,
        source_mapping_id=(
            PYSCF_ETF_NAC_MAPPING_ID_V233
            if use_etfs
            else PYSCF_NAC_MAPPING_ID_V233
        ),
        source_mapping_description=(
            PYSCF_NAC_EMPIRICAL_MAPPING_V232
            + (
                "; ETF translation removal enabled for a distinct diagnostic"
                if use_etfs
                else "; full many-electron overlap derivative"
            )
        ),
        use_etfs=use_etfs,
        mult_ediff=False,
    ).validate()


def analytic_nac_convention_v233(*, provider_name="analytic-provider"):
    provider_name = str(provider_name).strip()
    if not provider_name:
        raise ValueError("analytic provider name cannot be empty.")
    return DerivativeCouplingConventionV233(
        schema=NAC_CONVENTION_SCHEMA_V233,
        internal_definition=INTERNAL_NAC_DEFINITION_V233,
        source_backend=provider_name,
        source_backend_version="fixture",
        source_mapping_id=ANALYTIC_NAC_MAPPING_ID_V233,
        source_mapping_description=(
            "provider emits the framework internal d[i,j] convention directly"
        ),
        use_etfs=False,
        mult_ediff=False,
    ).validate()


def derivative_coupling_convention_from_dict_v233(payload):
    if not isinstance(payload, dict):
        raise TypeError("derivative-coupling convention payload must be a mapping.")
    expected = {
        "schema",
        "internal_definition",
        "source_backend",
        "source_backend_version",
        "source_mapping_id",
        "source_mapping_description",
        "use_etfs",
        "mult_ediff",
        "coordinate_unit",
    }
    if set(payload) != expected:
        raise ValueError("derivative-coupling convention field set mismatch.")
    return DerivativeCouplingConventionV233(**payload).validate()


@dataclass(frozen=True)
class LegacyReplayMigrationAttestationV233:
    legacy_dataset_fingerprint: str
    nac_disposition: str
    evidence: str

    def validate(self, *, expected_legacy_fingerprint, convention):
        convention = convention.validate()
        if self.legacy_dataset_fingerprint != str(expected_legacy_fingerprint):
            raise ValueError("legacy NAC attestation targets the wrong replay dataset.")
        if self.nac_disposition not in LEGACY_NAC_DISPOSITIONS_V233:
            raise ValueError("unknown legacy NAC disposition.")
        if not str(self.evidence).strip():
            raise ValueError("legacy NAC migration requires explicit evidence.")
        if self.nac_disposition in {"unknown", "requires_sign_correction"}:
            raise ValueError(
                "legacy replay NAC semantics are quarantined; automatic sign "
                "repair is forbidden."
            )
        if self.nac_disposition == "verified_corrected_v232":
            if convention.source_mapping_id not in {
                PYSCF_NAC_MAPPING_ID_V233,
                PYSCF_ETF_NAC_MAPPING_ID_V233,
            }:
                raise ValueError(
                    "corrected PySCF legacy data require the certified v0.23.2 "
                    "mapping identity."
                )
        if self.nac_disposition == "not_pyscf_derived":
            if convention.source_mapping_id != ANALYTIC_NAC_MAPPING_ID_V233:
                raise ValueError(
                    "non-PySCF legacy fixtures require an analytic convention."
                )
        return self

    def as_dict(self):
        return asdict(self)


def require_snapshot_nac_identity_v233(metadata, expected_convention):
    """Reject cache/snapshot metadata that omit or mismatch NAC identity."""
    if not isinstance(metadata, dict):
        raise TypeError("snapshot metadata must be a mapping.")
    expected = expected_convention.validate().fingerprint()
    observed = metadata.get("v233_nac_convention_fingerprint")
    if observed is None:
        raise ValueError(
            "snapshot lacks v0.23.3 NAC convention identity; legacy data are "
            "quarantined."
        )
    if observed != expected:
        raise ValueError("snapshot NAC convention identity mismatch.")
    return True
