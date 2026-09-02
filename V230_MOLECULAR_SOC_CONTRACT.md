# v0.23.0 molecular SOC contract

## Purpose

The v0.23.0 contract prevents a static SOC matrix, an analytic fixture, or an
untraceable file from being presented as a trajectory-ready molecular calculation.
The Gaussian engine still consumes the representation-neutral operators

$$
H=H_{\mathrm{sf}}+H_{\mathrm{SOC}},\qquad
K_a=K_{a,{\mathrm{sf}}}+K_{a,{\mathrm{SOC}}},\qquad
D_a=\langle\Phi|\partial_a\Phi\rangle.
$$

Admission concerns how those operators were obtained and validated, not a new
propagation equation.

## Capability tiers

`static_soc` requires a static SOC operator. It is suitable for point calculations and
must fail closed if asked to support moving-nuclear dynamics.

`trajectory_ready` is derived, not self-asserted. It requires all of:

1. static SOC;
2. spin-free physical derivatives;
3. SOC physical derivatives;
4. derivative connections;
5. cross-geometry overlaps.

Deterministic replay and analytic SOC derivatives are separate declarations. Declaring
analytic SOC derivatives without an SOC-derivative capability is invalid.

## Backend identity

Every provider freezes its backend name and version, source kind, electronic method,
basis, charge, electron count, SOC operator, scalar-relativistic method, derivative
method, active space, coordinate definition, state-tracking policy, and units. Internal
geometry, energy, and derivative units are exactly bohr, hartree, and hartree/bohr.

Allowed source kinds are:

- `validation_fixture`;
- `external_ab_initio_snapshot`;
- `live_ab_initio`.

Real sources additionally require a molecule name, ordered atom symbols, positive
isotope masses, an explicit reference geometry in bohr, and SHA-256 identities for the
calculation input and software environment. The electron-count parity must agree with
the SOC symmetry contract. This binds the molecular contract into the inherited
operator-provenance fingerprint.

## Evidence gates

A real source is not admitted by provenance alone. All five evidence families are
mandatory and independently evaluated:

| Gate | Required evidence |
|---|---|
| Independent reference | Named reference, measured error, and acceptance tolerance |
| Basis convergence | Ordered basis ladder and one change per adjacent pair |
| Method convergence | Ordered method ladder and one change per adjacent pair |
| Frame invariance | Translation and rotation residuals with a tolerance |
| State tracking | Minimum assigned overlap and assignment margin with thresholds |

Partial evidence declarations are invalid. A passing value in one family cannot
substitute for a missing value in another.

## Admission outcomes

The audit reports two distinct results:

- `protocol_passed`: the provider satisfies trajectory capabilities, convergence,
  parity/charge symmetry, component derivatives, cross-geometry differentials,
  provenance binding, and replay integrity;
- `real_backend_admitted`: every protocol check passes and the source is real,
  traceable, and supported by all five evidence families.

`require_molecular_soc_protocol_v230` is appropriate for deterministic validation.
`require_real_molecular_soc_backend_v230` is the only admission path for production
molecular SOC dynamics.

## Failure semantics

Admission fails on missing capabilities, false convergence, mixed parity or charge,
invalid time reversal or projectors, wrong component derivatives, wrong differential
connections, missing overlaps, changed provenance, corrupt replay data, incomplete
evidence, untraceable nuclei, or a fixture source. These are hard errors; there is no
warning-only production path.
