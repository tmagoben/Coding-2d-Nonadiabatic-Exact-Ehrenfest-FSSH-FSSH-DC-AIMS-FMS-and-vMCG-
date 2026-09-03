# v0.17 PySCF Error-Control Bridge

v0.17 is still validated on the analytic LVC model, but its error-control architecture
is directly relevant to a future molecular backend.

## 1. Why overlap-only locality is especially risky for ab initio dynamics

For molecular direct dynamics, matrix blocks can contain:

```text
electronic energy matrices
gradient-dependent terms
derivative couplings
local-diabatic transformations
moving-basis connection terms
```

A small nuclear overlap is useful screening information but cannot be assumed to bound
all of those terms uniformly.

The v0.17 S/H/T score is therefore a more appropriate architectural direction.

## 2. What would change for a molecular provider

The analytic release can evaluate $H_{ij}$ cheaply once the Gaussian pair data are
available.

For PySCF-backed dynamics, obtaining the electronic quantities needed for a candidate
edge can itself be expensive.

A molecular edge controller should therefore separate:

```text
cheap nuclear overlap/geometric screen
        ↓
electronic cache lookup
        ↓
provider-cost estimate
        ↓
only if necessary:
new electronic-structure calculation
        ↓
full edge importance
```

## 3. Audit hierarchy

A full dense molecular audit at every 20 steps would be prohibitively expensive.

v0.17 deliberately uses dense audits only because the analytic benchmark is cheap and
the goal is to establish the approximation hierarchy.

A molecular successor should replace the dense audit with a hierarchy such as:

```text
frequent:
    local omitted-edge score budget
    TDSE/projection residual diagnostics

less frequent:
    stratified omitted-edge samples
    selected dense graph neighborhoods

rare:
    full dense reference checkpoint on small validation systems
```

## 4. Gauge consistency

Any molecular $H_{ij}$ or derivative-coupling contribution used in the edge score must
be evaluated in consistent electronic frames.

The existing state-tracking and gauge-graph machinery from v0.6-v0.8 therefore remains
a prerequisite.

The nuclear sparse graph and electronic gauge graph remain conceptually distinct.

## 5. Provider-aware cost remains active

The v0.16 geometry-cache cost interface is retained.

A future candidate utility can therefore depend on both

$$
\text{predicted physical error reduction}
$$

and

$$
\text{electronic-structure cost}.
$$

This becomes particularly important before SOC, where electronic matrices may be
complex and the gauge structure more demanding.

## 6. What v0.17 still does not claim

v0.17 does not provide:

- a calibrated PySCF edge-error score;
- sampled molecular dense audits;
- asynchronous electronic-structure scheduling;
- sparse molecular TDSE defects;
- SOC matrix elements.

Those are intentionally deferred until the sparse approximation itself is better
controlled.
