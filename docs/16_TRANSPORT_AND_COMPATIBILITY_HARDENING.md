# Transport and compatibility hardening (v0.23.3)

v0.23.3 closes the last framework-level ambiguities identified after the real
spin-free PySCF validation:

1. Raw finite-manifold overlaps are contractions; electronic transport uses the
   separately certified unitary polar factor.
2. Physical overlap consistency and trajectory-quality thresholds are independent.
3. Replay format 2 binds raw overlaps, transports, singular values, policies, and
   derivative-coupling identity.
4. Legacy NAC data require explicit evidence. Unknown or wrong-sign records are
   quarantined without automatic repair.
5. Singlet/triplet and Kramers-doublet tracking operates on complete projector
   manifolds and is covariant under independent endpoint gauges.
6. Molecular-SOC matrix meaning is frozen down to prefactor, operator treatment,
   scalar relativity, exact state order, units, derivative semantics, and symmetry.
7. Provider identities bind every numerical convention used by replay, cache, or
   checkpoint paths.
8. Canonical byte identity and broader scientific compatibility are reported as
   separate runtime profiles.

These controls prepare the framework for one named method-specific molecular-SOC
source. They do not themselves validate or admit that source.
