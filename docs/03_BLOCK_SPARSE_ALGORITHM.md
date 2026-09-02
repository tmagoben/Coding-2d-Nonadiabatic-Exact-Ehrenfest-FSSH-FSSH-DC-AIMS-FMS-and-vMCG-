# Block-sparse molecular algorithm

## 1. Local basis

Each nuclear Gaussian \(g_i\) carries a complete \(s\)-dimensional local electronic frame,

$$
|\Xi_{i\alpha}\rangle=|g_i\rangle|\phi_{i\alpha}\rangle,
$$

so the total coefficient dimension is \(Ns\).

## 2. Cross-geometry overlap block

For electronic center frames \(\Phi_i\) and \(\Phi_j\),

$$
O_{ij}=\Phi_i^\dagger\Phi_j.
$$

With nuclear Gaussian overlap \(S_{ij}^{\rm nuc}\),

$$
\boxed{S_{ij}=S_{ij}^{\rm nuc}O_{ij}.}
$$

## 3. Pair-centroid Hamiltonian

At

$$
q_c=\frac12(q_i+q_j),
$$

transport both center frames into the centroid frame using unitary polar links \(U_{ci}\) and \(U_{cj}\). The electronic Hamiltonian contribution is

$$
H_{ij}^{e}=U_{ci}^\dagger H_e(q_c)U_{cj}.
$$

The block pair approximation is

$$
\boxed{
H_{ij}=T_{ij}^{\rm nuc}O_{ij}+S_{ij}^{\rm nuc}H_{ij}^{e}.
}
$$

This is the repository's discrete pair-centroid approximation; it is not claimed to be the complete production AIMS matrix element.

## 4. Moving-basis connection

Let

$$
\Gamma_j=\sum_aD_a(q_j)\dot q_{j,a}.
$$

Then the ordered seed is

$$
\boxed{
T_{ij}^{(0)}=\tau_{ij}^{\rm nuc}O_{ij}+S_{ij}^{\rm nuc}O_{ij}\Gamma_j.
}
$$

The reverse orientation is evaluated separately because \(T\) is not Hermitian.

## 5. Gauge-invariant sparse score

Define the normalized overlap magnitude

$$
s_{ij}=\frac{\|S_{ij}\|_F}{\sqrt s}.
$$

For the Hamiltonian,

$$
E_i=\frac{\|H_{ii}\|_F}{\sqrt s},
$$

$$
h_{ij}=\frac{\|H_{ij}\|_F/\sqrt s}{\max[\sqrt{E_iE_j},E_{\rm floor}]}.
$$

The nuclear moving-basis contribution is

$$
t_{ij}^{\rm nuc}=\Delta t\,
\frac{\sqrt{|\tau_{ij}|^2+|\tau_{ji}|^2}\,\|O_{ij}\|_F}{\sqrt{2s}}.
$$

For a truncated electronic subspace, use the singular values \(\sigma_k\) of \(O_{ij}\) to define

$$
m_{ij}=\left[\frac1s\sum_k(1-\sigma_k^2)\right]^{1/2}.
$$

The edge score is

$$
\boxed{
\eta_{ij}^2=(w_Ss_{ij})^2+(w_Hh_{ij})^2+(w_tt_{ij}^{\rm nuc})^2+(w_mm_{ij})^2.
}
$$

Each channel is invariant under independent unitary changes of the left and right electronic frames.

## 6. Hysteresis and accumulated omission

New edges require

$$
\eta_{ij}\ge\eta_{\rm enter},
$$

while existing edges persist until

$$
\eta_{ij}<\eta_{\rm exit}.
$$

For scored but omitted edges \(D\),

$$
B_D=\left(\sum_{e\in D}\eta_e^2\right)^{1/2}.
$$

If \(B_D\) exceeds the configured budget, the largest omitted edges are promoted until the budget is satisfied.

## 7. Sparse moving-basis propagation

The exact projected equation is

$$
iS\dot C=(H-iT)C.
$$

The metric identity is

$$
\dot S=T+T^\dagger.
$$

At each midpoint the implementation solves

$$
\boxed{
\left[S_m+\frac{\Delta t}{2}(iH_m+T_m)\right]C_{n+1}
=
\left[S_m-\frac{\Delta t}{2}(iH_m+T_m)\right]C_n.
}
$$

No inverse of \(S\) is formed.

## 8. Dynamic topology

The sparse graph is persistent by Gaussian uid and supports both edge entry and deletion as the Gaussian centers move. v0.21 includes a crossing/separation stress test specifically to ensure that the block propagator is not validated only on a frozen graph.
