# v0.24.2 program architecture

```text
PINNED PYSCF 2.13.1                    PROJECT SPIN ALGEBRA

7 independent OH calculations          spin-pure roots (E,S,Mref)
  q0 and q0 +/- {0.08,0.04,0.02}         complete |root,S,M_S> expansion
  ROHF -> equal 3-root SA-CASSCF          Wigner-Eckart / integer-2j CG
  common MOs + CI roots + density                    |
                 |                                   |
                 +----------------+------------------+
                                  v
                    DIRECT-JK BP-SOMF SNAPSHOT
  int1e_prinvxp                 int2e_p1vxp1 shell contractions
       |                       J - 3/2 K_L - 3/2 K_R
       +-----------------------+---------------------+
                               v
                 H_sf, H_SOC, H_total; time reversal
                 geometry + inputs + runtime + hashes
                               |
       explicit rank-five oracle at q0 only (validation, not production)
                               |
                               v
                 EXACT RESTRICTED-CASSCF OVERLAPS
                  O_root = <Psi(q0)|Psi(q+/-h)>
                               |
                    lift with delta_S delta_M
                               v
                 COMPLETE-MULTIPLET CONTRACTIONS O
                               |
                     O = U Sigma V^dagger
                               v
                 UNITARY POLAR TRANSPORT W = U V^dagger
                               |
             +-----------------+------------------+
             |                                    |
             v                                    v
 W H_sf(q+/-h) W^dagger               W H_SOC(q+/-h) W^dagger
             |                                    |
             +-----------------+------------------+
                               v
                    SEPARATE CENTERED DIFFERENCES
                       K_sf, K_SOC, K_total
                               |
                     3-step Richardson audit
                               |
              +----------------+-------------------+
              |                                    |
              v                                    v
 connected-geometry preview                  trajectory admission
 60/60 runtime gates                         FAIL CLOSED
 endpoint-bound evidence                     no full 3N tensor
 complete doublets                           no analytic derivative
                                              no physical continuous D
                                              no accuracy admission
```

## Trust chain

Each derivative record contains the center, minus, and plus snapshot fingerprints;
the raw cross-geometry contractions; both polar transports; four transported endpoint
component matrices; the three resulting derivatives; and residuals. The scan also
serializes compact receipts for all six endpoints. Validation recomputes endpoint
fingerprints, signed geometries, state order, runtime identity, state-average weights,
component derivatives, decomposition, Hermiticity, time reversal, and convergence.

The explicit rank-five SOMF tensor exists only in the center cross-check. No method
that exposes this differential preview is connected to the trajectory-provider
admission path.
