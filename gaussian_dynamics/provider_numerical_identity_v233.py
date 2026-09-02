"""Convention-complete provider identity for caches, replays, and checkpoints."""

from dataclasses import asdict, dataclass
import hashlib
import json

from .finite_manifold_transport_v233 import (
    OVERLAP_CONTRACT_ID_V233,
    TRANSPORT_CONTRACT_ID_V233,
    FiniteManifoldOverlapPolicyV233,
)
from .nac_compatibility_v233 import DerivativeCouplingConventionV233


PROVIDER_IDENTITY_SCHEMA_V233 = "gnd-provider-numerical-identity-v0.23.3"


def _canonical_bytes_v233(value):
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")


@dataclass(frozen=True)
class ProviderNumericalIdentityV233:
    schema: str
    provenance_fingerprint: str
    nac_convention_fingerprint: str
    overlap_contract: str
    transport_contract: str
    overlap_policy: dict
    replay_format_version: int

    def validate(self):
        if self.schema != PROVIDER_IDENTITY_SCHEMA_V233:
            raise ValueError("provider numerical identity schema mismatch.")
        for name in ("provenance_fingerprint", "nac_convention_fingerprint"):
            value = str(getattr(self, name))
            if len(value) != 64 or any(char not in "0123456789abcdef" for char in value):
                raise ValueError(f"{name} must be a lowercase SHA-256 digest.")
        if self.overlap_contract != OVERLAP_CONTRACT_ID_V233:
            raise ValueError("provider overlap contract is not v0.23.3.")
        if self.transport_contract != TRANSPORT_CONTRACT_ID_V233:
            raise ValueError("provider transport contract is not v0.23.3.")
        if int(self.replay_format_version) != 2:
            raise ValueError("provider replay format must be v0.23.3 format 2.")
        if not isinstance(self.overlap_policy, dict):
            raise TypeError("provider overlap policy must be a mapping.")
        expected = FiniteManifoldOverlapPolicyV233(**self.overlap_policy).validate()
        if expected.as_dict() != self.overlap_policy:
            raise ValueError("provider overlap policy is not canonical.")
        return self

    def as_dict(self):
        return asdict(self)

    def fingerprint(self):
        self.validate()
        return hashlib.sha256(_canonical_bytes_v233(self.as_dict())).hexdigest()


def build_provider_numerical_identity_v233(
    provenance,
    nac_convention,
    *,
    overlap_policy=FiniteManifoldOverlapPolicyV233(),
):
    provenance = provenance.validate()
    if type(nac_convention) is not DerivativeCouplingConventionV233:
        raise TypeError("nac_convention must be DerivativeCouplingConventionV233.")
    nac_convention = nac_convention.validate()
    overlap_policy = overlap_policy.validate()
    return ProviderNumericalIdentityV233(
        schema=PROVIDER_IDENTITY_SCHEMA_V233,
        provenance_fingerprint=provenance.fingerprint(),
        nac_convention_fingerprint=nac_convention.fingerprint(),
        overlap_contract=OVERLAP_CONTRACT_ID_V233,
        transport_contract=TRANSPORT_CONTRACT_ID_V233,
        overlap_policy=overlap_policy.as_dict(),
        replay_format_version=2,
    ).validate()


def require_provider_numerical_identity_v233(provider, expected_identity):
    if type(expected_identity) is not ProviderNumericalIdentityV233:
        raise TypeError("expected_identity must be ProviderNumericalIdentityV233.")
    expected_identity.validate()
    observed = getattr(provider, "numerical_identity_v233", None)
    if callable(observed):
        observed = observed()
    if type(observed) is not ProviderNumericalIdentityV233:
        raise ValueError(
            "provider lacks v0.23.3 numerical identity; legacy providers are "
            "quarantined from convention-bound caches and checkpoints."
        )
    observed.validate()
    if observed != expected_identity:
        raise ValueError("provider numerical identity differs from the trusted identity.")
    return observed


def run_convention_bound_dynamics_v233(
    provider,
    provenance,
    *,
    numerical_identity,
    **kwargs,
):
    """Run the inherited propagator with convention-complete checkpoint identity."""
    from .checkpoint_restart_v214 import run_self_consistent_block_dynamics_v214

    observed = require_provider_numerical_identity_v233(
        provider, numerical_identity
    )
    if observed.provenance_fingerprint != provenance.validate().fingerprint():
        raise ValueError(
            "provider numerical identity and supplied provenance disagree."
        )
    return run_self_consistent_block_dynamics_v214(
        provider,
        provenance,
        provider_numerical_fingerprint=observed.fingerprint(),
        **kwargs,
    )
