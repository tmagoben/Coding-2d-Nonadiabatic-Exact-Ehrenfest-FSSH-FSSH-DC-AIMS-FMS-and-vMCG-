# v0.17 Theory: Error-Controlled Sparse Gaussian Dynamics

v0.16 established a persistent overlap-locality graph and sparse projected dynamics.
Its remaining conceptual weakness was clear:

> A small Gaussian overlap is a rigorous statement about $S_{ij}$, but not by itself
> a rigorous upper bound on $H_{ij}$ or the moving-basis connection $T_{ij}$.

v0.17 therefore changes the locality decision itself.

The graph is now selected by an **exact local S/H/T importance score**, while overlap
geometry is retained only as a conservative search accelerator.

The release adds:

1. exact local S/H/T edge scoring;
2. score hysteresis;
3. a global local-importance L2 budget proxy;
4. periodic dense S/H/Snuc audits during propagation;
5. automatic one-sided threshold relaxation when an audit fails;
6. snapshot convergence campaigns in edge score and local score budget;
7. scaling tests using the actual S/H/T scorer rather than overlap locality alone.

## 1. Edge importance

For a candidate pair $(i,j)$ define

$$
s_{ij}=|S_{ij}|.
$$

Let

$$
h_{ij}
=
\frac{\|H_{ij}\|_F}
{\max\left[
\sqrt{\|H_{ii}\|_F\|H_{jj}\|_F},
E_{\rm floor}
\right]}.
$$

The time-connection contribution is

$$
t_{ij}
=
\Delta t
\sqrt{
|T_{ij}|^2+|T_{ji}|^2
}.
$$

The v0.17 dimensionless local score is

$$
\boxed{
\eta_{ij}
=
\sqrt{
(w_Ss_{ij})^2
+
(w_Hh_{ij})^2
+
(w_Tt_{ij})^2
}.
}
$$

The release weights are

$$
w_S=1,
\qquad
w_H=0.20,
\qquad
w_T=1.
$$

This is an explicit numerical importance metric. It is **not** claimed to equal a
rigorous operator-error bound.

## 2. Geometric search remains conservative

Before exact S/H/T scoring, v0.17 uses the v0.16 overlap upper bound

$$
|S_{ij}|
\le
\exp
\left[
-\frac12h_{ij}^{(q)}
\|q_i-q_j\|^2
\right].
$$

A safe KD-tree global radius is generated from the initial search floor

$$
\tau_{\rm search}=10^{-5}.
$$

Pairs inside that radius receive the tighter pair-specific overlap bound.

Only pairs surviving both geometric screens require exact pair algebra.

The geometric screen is therefore a **candidate-generation mechanism**. Final edge
membership is decided by the S/H/T score.

## 3. Score hysteresis

The deliberately aggressive initial release values are

$$
\eta_{\rm enter}=0.060,
\qquad
\eta_{\rm exit}=0.030.
$$

A new edge is admitted only above the enter score.

An existing edge persists until its score falls below the exit score.

## 4. Local importance budget

Individually small omitted edges can accumulate.

For all locally scored but omitted edges, define

$$
\boxed{
B_{\rm local}
=
\left(
\sum_{(i,j)\in D}
\eta_{ij}^2
\right)^{1/2}.
}
$$

The release imposes

$$
B_{\rm local}\le0.08.
$$

If the tentative dropped-edge set violates this budget, the largest omitted scores are
promoted back into the graph until the budget is satisfied.

This is a **global importance proxy**, not a proof that $\|\Delta H\|<0.08$.

## 5. Dense online audit

Every 20 steps the current sparse matrices are compared with a complete dense
v0.15-style pair build.

The audited quantities are

$$
\epsilon_S
=
\frac{\|S_s-S_d\|_F}{\|S_d\|_F},
$$

$$
\epsilon_H
=
\frac{\|H_s-H_d\|_F}{\|H_d\|_F},
$$

and

$$
\epsilon_{S_n}
=
\frac{\|S_{n,s}-S_{n,d}\|_F}
{\|S_{n,d}\|_F}.
$$

The release requires

$$
\boxed{
\epsilon_S,\epsilon_H,\epsilon_{S_n}
\le0.006.
}
$$

These audits are intentionally expensive. They are a correctness bridge for v0.17,
not the final large-$N$ strategy.

## 6. One-sided controller

If an audit fails, v0.17 conservatively relaxes both the S/H/T score thresholds and
the geometric search floor.

With relaxation factor $r=0.5$,

$$
\eta_{\rm enter}\leftarrow r\eta_{\rm enter},
$$

$$
\eta_{\rm exit}\leftarrow r\eta_{\rm exit},
$$

$$
\tau_{\rm search}\leftarrow r\tau_{\rm search}.
$$

Thresholds are **never tightened during the same run**.

## 7. Actual release controller behavior

The first audit intentionally fails:

```text
enter score: 0.060
S relative error: 0.01647878134803449
H relative error: 0.016940324970549453
```

The controller automatically changes to

```text
enter score: 0.03
exit score: 0.015
search overlap floor: 5e-06
```

and immediately passes:

```text
S relative error: 0.0026588687049987378
H relative error: 0.0020460833614795023
```

Every later scheduled audit also passes.

No unresolved audit occurs.

## 8. Final physical result

The final representation-consistent error is

$$
\boxed{
\|\rho_G-\rho_{\rm exact,projected}\|_F
=
0.0001336146005.
}
$$

Against the original exact target,

$$
\boxed{
\|\rho_G-\rho_{\rm target}\|_F
=
0.03333954068.
}
$$

The final v0.17 reduced density differs from v0.16 by only

$$
\boxed{
7.096e-15.
}
$$

In this benchmark the online controller therefore recovers the v0.16 accepted sparse
representation to floating-point precision.

## 9. Final matrix audit

At $t=0.6$,

$$
\epsilon_S=0.0051917427,
$$

$$
\epsilon_H=0.0039626323,
$$

and

$$
\epsilon_{S_n}=0.0051917427.
$$

All lie below the 0.006 budget.

## 10. Edge-threshold convergence

On the final snapshot, the score threshold is independently swept from

$$
0.12\rightarrow0.01.
$$

Both $S$ and $H$ errors decrease monotonically.

At the finest score threshold:

$$
\epsilon_S=0.00026621461,
$$

$$
\epsilon_H=0.00015704048.
$$

## 11. Local-budget convergence

Holding the nominal enter score fixed at $0.06$, the local omitted-score budget is
tightened until all candidate edges are restored.

Both audited $S$ and $H$ errors again decrease monotonically.

This demonstrates that the local importance budget is functioning as an error-control
coordinate rather than merely a diagnostic.

## 12. Adaptive basis growth remains intact

The physical TDSE-defect controller still adds one TBF at step 10.

The measured defect falls from

$$
0.032556131
$$

to

$$
0.030270054.
$$

The candidate enters with zero electronic amplitude, preserving wavefunction
continuity at the growth event.

## 13. Scientific interpretation

v0.17 establishes a layered control structure:

$$
\boxed{
\text{safe geometric search}
\rightarrow
\text{exact local S/H/T importance}
\rightarrow
\text{global local-score budget}
\rightarrow
\text{periodic dense matrix audit}
\rightarrow
\text{physical-grid TDSE defect}.
}
$$

That is a substantially stronger basis for later long-time convergence studies than an
overlap threshold alone.
